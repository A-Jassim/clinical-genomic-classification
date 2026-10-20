#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
for (module in c("data.R", "preprocessing.R", "modeling.R", "pipeline.R")) source(file.path(root, "R", module))
output <- commandArgs(trailingOnly = TRUE)
if (length(output) > 1) stop("Usage: Rscript scripts/run_analysis.R [output_directory]")
output <- if (length(output)) normalizePath(output, mustWork = FALSE) else file.path(root, "results")
run_pipeline(root, output)
