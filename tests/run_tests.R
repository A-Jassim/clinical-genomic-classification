#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
for (name in c("data.R", "preprocessing.R", "modeling.R", "pipeline.R")) source(file.path(root, "R", name))

stopifnot(identical(patient_id(c("TCGA-AA-0001-01A", "AA.0001.02B")), rep("TCGA-AA-0001", 2)))
stopifnot(auc_score(c(0, 0, 1, 1), c(0.1, 0.2, 0.8, 0.9)) == 1)
stopifnot(auc_score(c(0, 0, 1, 1), rep(0.5, 4)) == 0.5)
m <- classification_metrics(c(0, 0, 1, 1), c(0.1, 0.6, 0.7, 0.9), 0.5)
stopifnot(m$accuracy == 0.75, m$tp == 2, m$tn == 1, m$fp == 1, m$fn == 0)

training <- data.frame(age = c(20, 30, NA, 50), stage = c("T1", "T1", "T2", NA))
state <- fit_clinical_preprocessor(training, rare_count = 1)
assessment <- data.frame(age = c(NA, 999), stage = c("UNSEEN", NA))
encoded <- apply_clinical_preprocessor(assessment, state)
stopifnot(encoded[1, "age"] == 30, all(is.finite(encoded)), !any(grepl("UNSEEN", colnames(encoded))))
genes <- cbind(a = 1:4, b = 5:8)
gene_state <- fit_gene_preprocessor(genes)
new_genes <- matrix(c(NA, 10), nrow = 1, dimnames = list(NULL, c("a", "b")))
stopifnot(apply_gene_preprocessor(new_genes, gene_state)[1, "a"] == 2.5)

set.seed(80)
synthetic <- list(ids = sprintf("patient%03d", 1:90), y = rep(c(0L, 1L, 0L), 30),
                  clinical = data.frame(age = rnorm(90), stage = rep(c("T1", "T2"), 45)),
                  genes = matrix(rnorm(90 * 24), nrow = 90, dimnames = list(NULL, paste0("g", 1:24))))
synthetic$genes[, 1] <- synthetic$genes[, 1] + synthetic$y
split <- stratified_split(synthetic$y, 0.8, 2025)
stopifnot(length(intersect(split$train, split$test)) == 0, length(union(split$train, split$test)) == 90)
stopifnot(identical(split, stratified_split(synthetic$y, 0.8, 2025)))
config <- list(seed = 2025L, folds = 3L, top_genes = 8L, lambda_count = 8L)
first <- tune_models(synthetic, split$train, config)
# Changing every held-out feature and outcome must not affect any tuning result.
changed <- synthetic
changed$y[split$test] <- 1L - changed$y[split$test]
changed$genes[split$test, ] <- 1e6
changed$clinical$age[split$test] <- -1e6
changed$clinical$stage[split$test] <- "NEW"
second <- tune_models(changed, split$train, config)
stopifnot(identical(first, second))
stopifnot(all(!synthetic$ids[split$test] %in% first$oof$patient_id))
for (i in seq_len(nrow(first$selected))) {
  chosen <- first$selected[i, ]
  oof <- first$oof[first$oof$model == chosen$model, ]
  stopifnot(chosen$threshold == select_threshold(oof$actual, oof$probability))
}
message("All tests passed: identity normalization, metrics, training-only preprocessing, deterministic splitting, and held-out-data isolation.")
