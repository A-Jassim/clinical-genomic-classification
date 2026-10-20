#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
for (module in c("data.R", "preprocessing.R", "modeling.R", "pipeline.R")) source(file.path(root, "R", module))
result <- readRDS(file.path(root, "results", "analysis.rds"))
partition <- result$partitions
stopifnot(!anyDuplicated(partition$patient_id))
train <- partition$patient_id[partition$partition == "train"]
test <- partition$patient_id[partition$partition == "test"]
stopifnot(length(intersect(train, test)) == 0)
stopifnot(setequal(result$out_of_fold$patient_id, train))
stopifnot(setequal(result$predictions$patient_id, test))
stopifnot(all(is.na(partition$fold[partition$partition == "test"])))
stopifnot(result$selected_model == result$validation$model[which.max(result$validation$validation_auc)])
for (name in result$metrics$model) {
  validation <- result$out_of_fold[result$out_of_fold$model == name, ]
  summary <- result$metrics[result$metrics$model == name, ]
  threshold <- select_threshold(validation$actual, validation$probability)
  stopifnot(isTRUE(all.equal(threshold, summary$threshold)))
  predictions <- result$predictions[result$predictions$model == name, ]
  stopifnot(nrow(predictions) == length(test), !anyDuplicated(predictions$patient_id))
  recomputed <- classification_metrics(predictions$actual, predictions$probability, threshold)
  stopifnot(isTRUE(all.equal(as.numeric(summary[, names(recomputed)]), as.numeric(recomputed))))
}
message("Export integrity, patient separation, training-derived thresholds, and metrics verified.")
if ("--repeat" %in% commandArgs(trailingOnly = TRUE)) {
  repeated <- run_pipeline(root, file.path(root, "results", "reproducibility_check"), result$config)
  for (field in c("metrics", "predictions", "partitions", "validation", "validation_grid", "out_of_fold", "fold_screening", "selected_genes")) {
    if (!identical(result[[field]], repeated[[field]])) stop("Reproduction mismatch: ", field)
  }
  message("Full rerun reproduced metrics, predictions, partitions, validation, and selected genes exactly.")
}
