#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
if (!file.exists(file.path(root, "results", "analysis.rds"))) stop("Run scripts/run_analysis.R first.")
if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available()) stop("R Markdown and Pandoc are required for rendering.")
dir.create(file.path(root, "results", "reports"), recursive = TRUE, showWarnings = FALSE)
for (name in c("comprehensive_analysis_and_model_comparison.Rmd", "weighted_genomic_clinical_modeling.Rmd")) {
  rmarkdown::render(file.path(root, "notebooks", name), knit_root_dir = root,
                     output_dir = file.path(root, "results", "reports"), envir = new.env(parent = globalenv()), quiet = TRUE)
}
