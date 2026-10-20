# Load raw objects into an isolated environment and align by patient identity.
patient_id <- function(x) {
  x <- toupper(gsub(".", "-", trimws(as.character(x)), fixed = TRUE))
  x <- ifelse(grepl("^TCGA-", x), x, paste0("TCGA-", x))
  valid <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}($|-)", x)
  if (any(!valid)) stop("Unrecognized TCGA patient identifiers: ", paste(head(x[!valid]), collapse = ", "))
  substr(x, 1, 12)
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | tolower(x) %in% c("", "na", "nan", "unknown", "not reported", "missing", "[not available]", "[not applicable]")
  x[missing] <- NA_character_
  x
}

clinical_features <- function(clinical) {
  n <- nrow(clinical)
  get <- function(name) if (name %in% names(clinical)) clinical[[name]] else rep(NA_character_, n)
  diagnosis <- clean_text(get("primary_diagnosis"))
  grouped <- rep("Other", n)
  patterns <- c(Ductal = "Infiltrating duct|Ductal", Lobular = "Lobular", Mucinous = "Mucinous|Colloid", Medullary = "Medullary", Papillary = "Papillary", Metaplastic = "Metaplastic")
  for (name in names(patterns)) grouped[!is.na(diagnosis) & grepl(patterns[[name]], diagnosis, ignore.case = TRUE)] <- name
  grouped[is.na(diagnosis)] <- NA_character_
  stage <- toupper(clean_text(get("ajcc_pathologic_t")))
  for (k in 1:4) stage[!is.na(stage) & grepl(paste0("^T", k), stage)] <- paste0("T", k)
  sites <- get("sites_of_involvement")
  if (!is.list(sites)) sites <- as.list(sites)
  sites <- lapply(sites, function(x) { x <- clean_text(unlist(x)); x[!is.na(x)] })
  absent <- lengths(sites) == 0
  count <- function(pattern) vapply(sites, function(x) sum(grepl(pattern, x, ignore.case = TRUE)), numeric(1))
  out <- data.frame(
    age_at_index = suppressWarnings(as.numeric(as.character(get("age_at_index")))),
    initial_weight = suppressWarnings(as.numeric(as.character(get("initial_weight")))),
    laterality = clean_text(get("laterality")),
    ajcc_pathologic_t = stage,
    prior_treatment = clean_text(get("prior_treatment")),
    diagnosis_grouped = grouped,
    sites_data_missing = as.numeric(absent),
    nb_sites = as.numeric(lengths(sites)),
    bilateral = as.numeric(count("Left") > 0 & count("Right") > 0),
    nb_quadrants = count("Upper Outer|Upper Inner|Lower Outer|Lower Inner"),
    has_metastasis = as.numeric(count("Bone|Lung|Liver|Lymph node|Skin|Ovary") > 0),
    stringsAsFactors = FALSE
  )
  # Missing site information is unknown, not evidence of absence.
  for (name in c("nb_sites", "bilateral", "nb_quadrants", "has_metastasis")) out[absent, name] <- NA_real_
  # Follow-up duration and treatment response are deliberately excluded.
  out
}

load_project_data <- function(path) {
  if (!requireNamespace("S4Vectors", quietly = TRUE)) stop("Install S4Vectors first: Rscript scripts/setup.R")
  env <- new.env(parent = baseenv())
  load(path, envir = env)
  if (!all(c("GeneX", "clinical_data") %in% ls(env))) stop("Raw file must contain GeneX and clinical_data.")
  genes <- as.matrix(env$GeneX)
  storage.mode(genes) <- "double"
  clinical <- as.data.frame(env$clinical_data)
  raw_clinical_rows <- nrow(clinical)
  if (is.null(rownames(genes)) || is.null(colnames(genes))) stop("GeneX needs row and column identifiers.")
  if (anyDuplicated(colnames(genes))) stop("Duplicate gene names.")
  clinical_ids <- if ("bcr_patient_barcode" %in% names(clinical)) clinical$bcr_patient_barcode else rownames(clinical)
  gene_ids <- patient_id(rownames(genes))
  clinical_ids <- patient_id(clinical_ids)
  raw_status <- tolower(clean_text(clinical$vital_status))
  y <- rep(NA_integer_, length(raw_status))
  y[raw_status %in% c("alive", "living", "0")] <- 0L
  y[raw_status %in% c("dead", "deceased", "1")] <- 1L
  unknown <- !is.na(raw_status) & is.na(y)
  if (any(unknown)) stop("Unrecognized vital status: ", paste(unique(raw_status[unknown]), collapse = ", "))
  # Refuse contradictory records rather than arbitrarily selecting a target.
  for (id in unique(clinical_ids[duplicated(clinical_ids)])) {
    values <- unique(na.omit(y[clinical_ids == id]))
    if (length(values) > 1) stop("Conflicting vital status for patient ", id)
  }
  # Prefer a labeled record; ties retain original input order deterministically.
  order_clinical <- order(is.na(y), seq_along(y))
  keep_clinical <- order_clinical[!duplicated(clinical_ids[order_clinical])]
  clinical <- clinical[keep_clinical, , drop = FALSE]
  y <- y[keep_clinical]
  clinical_ids <- clinical_ids[keep_clinical]
  # At most one gene-expression sample per patient, independent of the target.
  keep_genes <- !duplicated(gene_ids)
  dropped_samples <- sum(!keep_genes)
  genes <- genes[keep_genes, , drop = FALSE]
  gene_ids <- gene_ids[keep_genes]
  ids <- sort(intersect(gene_ids, clinical_ids[!is.na(y)]))
  if (length(ids) < 50) stop("Fewer than 50 aligned, labeled patients.")
  gi <- match(ids, gene_ids)
  ci <- match(ids, clinical_ids)
  rownames(genes) <- gene_ids
  features <- clinical_features(clinical[ci, , drop = FALSE])
  rownames(features) <- ids
  list(ids = ids, genes = genes[gi, , drop = FALSE], clinical = features, y = y[ci],
       audit = list(raw_clinical_rows = raw_clinical_rows,
                    aligned_patients = length(ids), dropped_extra_gene_samples = dropped_samples,
                    genes = ncol(genes), alive = sum(y[ci] == 0), dead = sum(y[ci] == 1)))
}

stratified_split <- function(y, fraction, seed) {
  set.seed(seed)
  train <- sort(unlist(lapply(0:1, function(k) {
    indices <- which(y == k)
    if (length(indices) < 4) stop("Not enough patients in each class.")
    sample(indices, max(2, min(length(indices) - 2, floor(length(indices) * fraction))))
  })))
  list(train = train, test = setdiff(seq_along(y), train))
}

stratified_folds <- function(y, k, seed) {
  if (min(table(y)) < k) stop("Too few minority-class patients for requested folds.")
  set.seed(seed)
  folds <- integer(length(y))
  for (value in 0:1) {
    idx <- which(y == value)
    folds[sample(idx)] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}
