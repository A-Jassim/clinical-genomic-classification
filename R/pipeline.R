default_config <- function() list(seed = 2025L, train_fraction = 0.8, folds = 5L, top_genes = 100L, lambda_count = 40L)

run_pipeline <- function(root, output = file.path(root, "results"), config = default_config()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install dependencies first: Rscript scripts/setup.R")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(config$seed)
  # Keep explicit lambda grids intact for reproducible comparison across folds.
  old_control <- glmnet::glmnet.control()
  on.exit(do.call(glmnet::glmnet.control, old_control), add = TRUE)
  glmnet::glmnet.control(fdev = 0, devmax = 1)
  input <- file.path(root, "mrr_bio.Rdata")
  data <- load_project_data(input)
  split <- stratified_split(data$y, config$train_fraction, config$seed)
  stopifnot(length(intersect(data$ids[split$train], data$ids[split$test])) == 0)
  message("Aligned ", length(data$ids), " unique patients: ", length(split$train), " training / ", length(split$test), " test.")
  tuning <- tune_models(data, split$train, config)
  # Freeze model choice before looking at any test predictions or metrics.
  winner <- tuning$selected$model[which.max(tuning$selected$validation_auc)]
  combined_rows <- which(tuning$selected$features == "combined")
  combined_winner <- tuning$selected$model[combined_rows[which.max(tuning$selected$validation_auc[combined_rows])]]
  prepared <- prepare_partitions(data, split$train, split$test, config$top_genes)
  summaries <- predictions <- vector("list", nrow(tuning$selected))
  saved_models <- list()
  for (i in seq_len(nrow(tuning$selected))) {
    candidate <- tuning$selected[i, ]
    path <- fit_path(prepared$train[[candidate$features]], data$y[split$train], candidate,
                     lambda_grid(candidate$alpha, config$lambda_count))
    probability <- as.numeric(predict(path, newx = prepared$assessment[[candidate$features]], s = candidate$lambda, type = "response"))
    summaries[[i]] <- cbind(candidate, selected_by_validation = candidate$model == winner,
                             classification_metrics(data$y[split$test], probability, candidate$threshold))
    predictions[[i]] <- data.frame(patient_id = data$ids[split$test], model = candidate$model,
                                   actual = data$y[split$test], probability = probability,
                                   threshold = candidate$threshold, predicted = as.integer(probability >= candidate$threshold))
    if (candidate$model %in% c(winner, combined_winner)) {
      saved_models[[candidate$model]] <- list(fit = path, candidate = candidate, preprocessing = prepared$state)
    }
  }
  metrics <- do.call(rbind, summaries)
  predictions <- do.call(rbind, predictions)
  assignment <- data.frame(patient_id = data$ids, actual = data$y,
                            partition = ifelse(seq_along(data$ids) %in% split$train, "train", "test"), fold = NA_integer_)
  assignment$fold[split$train] <- tuning$folds
  # Descriptive summaries use training observations only.
  clinical_train <- data$clinical[split$train, , drop = FALSE]
  missingness <- data.frame(feature = names(clinical_train), missing_fraction = vapply(clinical_train, function(x) mean(is.na(x)), numeric(1)))
  result <- list(config = config, input_md5 = unname(tools::md5sum(input)), data_audit = data$audit,
                 selected_model = winner, selected_combined_model = combined_winner,
                 metrics = metrics, predictions = predictions, partitions = assignment,
                 validation = tuning$selected, validation_grid = tuning$grid,
                 out_of_fold = tuning$oof, fold_screening = tuning$screening,
                 selected_genes = prepared$state$selected_genes, missingness = missingness,
                 training_class_counts = table(factor(data$y[split$train], levels = 0:1, labels = c("Alive", "Dead"))),
                 package_versions = vapply(c("glmnet", "Matrix", "S4Vectors", "jsonlite", "rmarkdown", "knitr"), function(p) as.character(utils::packageVersion(p)), character(1)))
  saveRDS(result, file.path(output, "analysis.rds"), version = 3)
  saveRDS(saved_models, file.path(output, "fitted_models.rds"), version = 3)
  write.csv(metrics, file.path(output, "test_metrics.csv"), row.names = FALSE)
  write.csv(predictions, file.path(output, "test_predictions.csv"), row.names = FALSE)
  write.csv(assignment, file.path(output, "patient_partitions.csv"), row.names = FALSE)
  write.csv(tuning$selected, file.path(output, "validation_selection.csv"), row.names = FALSE)
  write.csv(tuning$grid, file.path(output, "validation_grid.csv"), row.names = FALSE)
  write.csv(tuning$oof, file.path(output, "out_of_fold_predictions.csv"), row.names = FALSE)
  write.csv(tuning$screening, file.path(output, "fold_gene_screening.csv"), row.names = FALSE)
  write.csv(data.frame(gene = result$selected_genes), file.path(output, "selected_genes.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  write_results_summary(result, output)
  jsonlite::write_json(list(config = config, input_md5 = result$input_md5,
                           data_audit = result$data_audit, selected_model = winner,
                           selected_combined_model = combined_winner,
                           package_versions = as.list(result$package_versions)),
                       file.path(output, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  message("Frozen training-selected model: ", winner)
  print(metrics[, c("model", "validation_auc", "auc", "accuracy", "sensitivity", "specificity")], row.names = FALSE)
  invisible(result)
}

write_results_summary <- function(result, output) {
  chosen <- result$metrics[result$metrics$model == result$selected_model, ]
  test <- result$predictions[result$predictions$model == result$selected_model, ]
  baseline <- mean(test$actual == 0)
  rows <- apply(result$metrics, 1, function(row) sprintf(
    "| %s | %.3f | %.3f | %.1f%% | %.1f%% | %.1f%% |",
    row[["model"]], as.numeric(row[["validation_auc"]]), as.numeric(row[["auc"]]),
    100 * as.numeric(row[["accuracy"]]), 100 * as.numeric(row[["sensitivity"]]), 100 * as.numeric(row[["specificity"]])))
  text <- c(
    "# Recomputed results", "",
    "Generated by the corrected shared pipeline. These results replace neither the text nor the figures of the original academic PDF reports.", "",
    sprintf("- Unique aligned patients: **%d** (%d Alive, %d Dead).", result$data_audit$aligned_patients, result$data_audit$alive, result$data_audit$dead),
    sprintf("- Train/test: **%d / %d**, with a patient-level stratified split and seed %d.", sum(result$partitions$partition == "train"), nrow(test), result$config$seed),
    sprintf("- Test composition: **%d Alive / %d Dead**.", sum(test$actual == 0), sum(test$actual == 1)),
    sprintf("- Training-selected model: **%s**; lambda **%.6g**; frozen threshold **%.4f**.", result$selected_model, chosen$lambda, chosen$threshold), "",
    "## Held-out test metrics", "",
    "| Model | Tuning AUC | Test AUC | Accuracy | Sensitivity | Specificity |",
    "| --- | ---: | ---: | ---: | ---: | ---: |", rows, "",
    sprintf("The training-selected model has %d true positives, %d true negatives, %d false positives, and %d false negatives on the test set.", chosen$tp, chosen$tn, chosen$fp, chosen$fn), "",
    sprintf("**Class-imbalance reference:** predicting Alive for every test patient would achieve %.1f%% accuracy, 0%% sensitivity for Dead, and AUC 0.5. The chosen model does not outperform this trivial baseline on accuracy; its AUC and minority-class detection answer different questions. Accuracy alone is insufficient here.", 100 * baseline), "",
    "The selected model is determined by training validation AUC, not by test ranking. Validation AUC is a tuning score. The single held-out split is not an external validation study, and the gene panel is exploratory.", "",
    "## Reproduction", "",
    "Run `Rscript scripts/run_analysis.R`, then `Rscript scripts/verify_results.R --repeat` to verify output integrity and compare a second complete run.", "",
    sprintf("Input MD5: `%s`. R and package versions are recorded in `renv.lock` and `results/sessionInfo.txt`.", result$input_md5)
  )
  writeLines(text, file.path(output, "RESULTS.md"), useBytes = TRUE)
  grDevices::png(file.path(output, "roc_comparison.png"), width = 1400, height = 900, res = 150)
  tryCatch(plot_roc_comparison(result), finally = grDevices::dev.off())
}

plot_roc_comparison <- function(result, models = result$metrics$model) {
  colors <- grDevices::hcl.colors(length(models), "Dark 3")
  plot(c(0, 1), c(0, 1), type = "n", xlab = "False positive rate", ylab = "Sensitivity",
       main = "Held-out test ROC curves", xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  abline(0, 1, col = "grey70", lty = 2)
  for (i in seq_along(models)) {
    p <- result$predictions[result$predictions$model == models[i], ]
    roc <- roc_points(p$actual, p$probability)
    lines(1 - roc$specificity, roc$sensitivity, col = colors[i], lwd = 2)
  }
  auc <- result$metrics$auc[match(models, result$metrics$model)]
  legend("bottomright", legend = sprintf("%s (AUC %.3f)", models, auc), col = colors, lwd = 2, cex = 0.65, bty = "n")
}
