# Verification record

The corrected workflow was executed on the original dataset in Debian using R 4.2.2.

- The training/test patient identities are disjoint.
- Out-of-fold predictions cover training patients only; test predictions cover test patients only.
- Thresholds were recomputed from training out-of-fold predictions and matched the exported values.
- Metrics were independently recomputed from exported test predictions.
- Synthetic tests changed every held-out outcome and feature without changing training-side tuning results.
- Tests cover unseen clinical categories, training-median imputation (including genes), deterministic splits, patient-ID normalization, and tied-score AUC.
- A complete rerun in the final locked environment reproduced all eight exported CSVs byte for byte, including predictions, metrics, patient assignments, and gene selection.
- Both HTML notebooks compiled successfully with the recorded dependencies.
- `renv::restore()` found the library synchronized; `renv::status()` reported a consistent project.
- A static `codetools::checkUsage()` check reported no issues in the shared functions.

These checks establish reproducibility on the tested dataset and environment, not external clinical validity or guaranteed performance on other cohorts.
