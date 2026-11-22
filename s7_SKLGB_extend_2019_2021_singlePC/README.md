# Extended 2019–2021 test-year evaluation (one PC)

Retains the fixed 2,000-cell stratified sample and evaluates the original nine forecasting configurations using new retrospective test origins in 2019, 2020 and 2021. For test year T, the intended rolling design uses training through T-2, validation at T-1, and test at T.

## Inputs

- `SKLGB_hyperparameter_sample2000.mat`, included, is the same fixed sample used in other sensitivity analyses.
- The full original `GLDAS_common_valid_mask_1deg_200401_202412.mat` must be available under the source project root; it is not included.
- MATLAB System Identification, Econometrics, Statistics and Machine Learning, Deep Learning and optional Parallel Computing toolboxes.

## Instructions

1. Edit `cfg.projectRoot` (and parallel settings, if needed) in `SKLGB_ext_config.m`.
2. In MATLAB, add this folder to the path and execute:

```matlab
SKLGB_ext_run('preflight');
SKLGB_ext_smoke('AR');
SKLGB_ext_smoke('ARIMA');
SKLGB_ext_smoke('SARIMA');
SKLGB_ext_run('train');
SKLGB_ext_run('analyze');
```

The training runner executes the nine configured models sequentially and can skip completed compatible outputs. Linear fits use the bundled original-source-derived `SKLGB_ext_*_core.m` functions, with SHA-256 provenance checks for existing results. The folder `Original_linear_sources/` contained archival copies of the original global scripts rather than runtime dependencies; these duplicates were removed.

Expected summary files are placed under `<projectRoot>/s3_darft/s1_analysis/s7_extended_2019_2021_2000grids_singlePC/annual_analysis/`. This includes early-vs-late descriptive comparisons and exploratory paired spatial-block bootstrap intervals. The bootstrap does not imply independent replication of years.
