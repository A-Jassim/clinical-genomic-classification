#!/usr/bin/env Rscript
# Maintainer utility: record an intentionally updated, tested environment.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
packages <- c("glmnet", "jsonlite", "S4Vectors", "rmarkdown", "knitr", "BiocManager", "renv")
if (!requireNamespace("renv", quietly = TRUE)) stop("Install renv before updating the environment.")
renv::load(project = root)
renv::snapshot(project = root, packages = packages, prompt = FALSE)
