# Each learned transformation is fitted on a training partition only.
fit_clinical_preprocessor <- function(x, rare_count = 5L) {
  specification <- lapply(x, function(column) {
    if (is.numeric(column)) {
      finite <- column[is.finite(column)]
      return(list(type = "numeric", median = if (length(finite)) median(finite) else 0))
    }
    column <- clean_text(column)
    counts <- table(column)
    list(type = "categorical", known = sort(names(counts)[counts >= rare_count]))
  })
  state <- list(specification = specification, columns = NULL)
  matrix <- apply_clinical_preprocessor(x, state)
  state$columns <- colnames(matrix)[apply(matrix, 2, function(column) length(unique(column)) > 1)]
  if (length(state$columns) < 2) stop("Fewer than two varying clinical predictors.")
  state
}

apply_clinical_preprocessor <- function(x, state) {
  pieces <- lapply(names(state$specification), function(name) {
    spec <- state$specification[[name]]
    column <- x[[name]]
    if (spec$type == "numeric") {
      missing <- !is.finite(column)
      column[missing] <- spec$median
      matrix <- cbind(column, as.numeric(missing))
      colnames(matrix) <- paste0(name, c("", "__missing"))
      return(matrix)
    }
    column <- clean_text(column)
    column[!is.na(column) & !column %in% spec$known] <- "__OTHER__"
    column[is.na(column)] <- "__MISSING__"
    levels <- c(spec$known, "__MISSING__", "__OTHER__")
    matrix <- outer(column, levels, "==") * 1
    colnames(matrix) <- paste0(name, "=", levels)
    matrix
  })
  matrix <- do.call(cbind, pieces)
  if (!is.null(state$columns)) matrix <- matrix[, state$columns, drop = FALSE]
  if (any(!is.finite(matrix))) stop("Nonfinite clinical design matrix.")
  matrix
}

fit_gene_preprocessor <- function(x) {
  medians <- apply(x, 2, function(column) {
    values <- column[is.finite(column)]
    if (length(values)) median(values) else 0
  })
  for (j in which(colSums(!is.finite(x)) > 0)) {
    x[!is.finite(x[, j]), j] <- medians[j]
  }
  variances <- colSums((sweep(x, 2, colMeans(x)))^2)
  genes <- colnames(x)[variances > 1e-12]
  if (length(genes) < 2) stop("Fewer than two varying genes.")
  list(medians = medians[genes], genes = genes)
}

apply_gene_preprocessor <- function(x, state) {
  x <- x[, state$genes, drop = FALSE]
  for (j in which(colSums(!is.finite(x)) > 0)) x[!is.finite(x[, j]), j] <- state$medians[j]
  x
}

screen_genes <- function(x, y, top_n) {
  # Welch statistics; selection is recomputed inside every validation fold.
  a <- x[y == 0, , drop = FALSE]
  b <- x[y == 1, , drop = FALSE]
  variance <- function(z) colSums(sweep(z, 2, colMeans(z))^2) / (nrow(z) - 1)
  va <- variance(a) / nrow(a)
  vb <- variance(b) / nrow(b)
  statistic <- abs(colMeans(a) - colMeans(b)) / sqrt(va + vb)
  df <- (va + vb)^2 / (va^2 / (nrow(a) - 1) + vb^2 / (nrow(b) - 1))
  p <- 2 * pt(-statistic, df = df)
  p[!is.finite(p)] <- 1
  colnames(x)[order(p, colnames(x))[seq_len(min(top_n, ncol(x)))]]
}

prepare_partitions <- function(data, train, assessment, top_genes) {
  clinical_state <- fit_clinical_preprocessor(data$clinical[train, , drop = FALSE])
  gene_state <- fit_gene_preprocessor(data$genes[train, , drop = FALSE])
  c_train <- apply_clinical_preprocessor(data$clinical[train, , drop = FALSE], clinical_state)
  c_assess <- apply_clinical_preprocessor(data$clinical[assessment, , drop = FALSE], clinical_state)
  g_train <- apply_gene_preprocessor(data$genes[train, , drop = FALSE], gene_state)
  g_assess <- apply_gene_preprocessor(data$genes[assessment, , drop = FALSE], gene_state)
  selected <- screen_genes(g_train, data$y[train], top_genes)
  colnames(c_train) <- paste0("clinical::", colnames(c_train))
  colnames(c_assess) <- colnames(c_train)
  colnames(g_train) <- paste0("gene::", colnames(g_train))
  colnames(g_assess) <- colnames(g_train)
  selected_columns <- paste0("gene::", selected)
  list(train = list(clinical = c_train, genomic = g_train,
                    combined = cbind(c_train, g_train[, selected_columns, drop = FALSE])),
       assessment = list(clinical = c_assess, genomic = g_assess,
                         combined = cbind(c_assess, g_assess[, selected_columns, drop = FALSE])),
       state = list(clinical = clinical_state, genomic = gene_state, selected_genes = selected))
}
