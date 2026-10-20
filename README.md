# Clinical and Genomic Classification with Regularized Regression

**ENSIIE · MRR academic project · R**  


Compare clinical features, 5,000 gene-expression measurements, and their combination to classify observed patient vital status (`Alive` / `Dead`) in a TCGA-derived dataset. This repository provides a single executable analysis and two complementary reports built from its outputs.

## Reproduce the study

Run these commands from the repository root with R available:

```sh
Rscript scripts/setup.R
Rscript tests/run_tests.R
Rscript scripts/run_analysis.R
Rscript scripts/verify_results.R
Rscript scripts/render_reports.R
```

R Markdown reports are written to `results/reports/`. The analysis also exports test metrics, predictions, patient partitions, validation folds, selected genes, fitted model objects, and session information. For a full deterministic rerun check:

```sh
Rscript scripts/verify_results.R --repeat
```

Use the R version recorded in `renv.lock` for the closest reproduction. `renv::restore()` restores recorded R package versions; it does not install R itself or operating-system tools. Pandoc is needed only for HTML report rendering. On Debian, install R, a C++ toolchain and Pandoc before restoring packages:

```sh
sudo apt-get update
sudo apt-get install r-base r-base-dev pandoc
```

`scripts/setup.R` restores the package versions recorded in `renv.lock`. The raw clinical object is a Bioconductor `S4Vectors` object, so that dependency is required even though preprocessing uses base R after conversion.

## One pipeline, two reports

| Entry point | Purpose |
| --- | --- |
| `scripts/run_analysis.R` | Fit and evaluate every candidate using one shared protocol |
| [Complete analysis notebook](notebooks/comprehensive_analysis_and_model_comparison.Rmd) | Training-data exploration, protocol, validation, and test comparison |
| [Weighted modeling notebook](notebooks/weighted_genomic_clinical_modeling.Rmd) | Examine the combined candidate selected using training validation |

The notebooks read `results/analysis.rds`; they do not duplicate fitting logic or create alternative splits. They are rendered with `scripts/render_reports.R`, which supplies the correct computation directory.

## Validation protocol

1. Load only `mrr_bio.Rdata`, align records by normalized TCGA patient identity, exclude unknown outcomes, and retain one expression sample per patient deterministically. Conflicting duplicate outcomes cause a clear error.
2. Reserve 20% of unique patients as a stratified test set using seed 2025.
3. Use five stratified folds within the training partition. Fit numerical imputation, categorical encoding, rare-level handling, gene variance filtering, and gene screening independently within each fold.
4. Compare clinical Ridge/Lasso/Elastic Net, genomic Ridge/Lasso/Elastic Net, weighted genomic Lasso, and weighted combined Ridge/Lasso/Elastic Net. Combined models use the top 100 genes by training-fold Welch-test p-value; this is screening, not a relaxed-Lasso implementation.
5. Choose lambda and model family by pooled out-of-fold AUC. Choose each threshold by out-of-fold Youden index. Both are frozen before final test evaluation.
6. Refit each candidate on the full training partition and evaluate on the same held-out patients. The model with the best test score does not replace the training-selected model.

Numeric medians, categorical levels, class weights, and glmnet standardization are all fitted without test observations. The target is never used as an imputation predictor. Follow-up duration, treatment response, and patient identifiers are excluded from clinical predictors.

## Structure

```text
├── R/
│   ├── data.R                 # Identity alignment, feature extraction, splits
│   ├── preprocessing.R        # Training-only transformations and gene screening
│   ├── modeling.R             # Candidate models, validation, thresholds, metrics
│   └── pipeline.R             # Shared fitting workflow and report helpers
├── scripts/                   # Setup, execution, verification, and rendering
├── tests/run_tests.R          # Leakage-isolation and preprocessing checks
├── notebooks/                # Two reports using shared pipeline outputs
├── reports/                  # Original academic PDF reports and figures
├── results/                   # Generated analysis outputs and reports
├── renv.lock                  # Recorded R package versions
└── mrr_bio.Rdata              # Original input; never overwritten
```

The old `data_clean.csv` is no longer needed. Clinical features are rebuilt from the raw input for every run. The earlier notebooks were backed up outside the repository before being replaced.

## Data and attribution

The original project attributes the supplied data to **The Cancer Genome Atlas (TCGA)** and cites *Comprehensive molecular portraits of human breast tumours* (The Cancer Genome Atlas Network, 2012), DOI `10.1038/nature11412`. The raw dataset and original reports are preserved. Exact data extraction and redistribution terms were not documented in the supplied folder.