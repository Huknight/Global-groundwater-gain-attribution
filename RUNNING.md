# Running the code: overview and dependencies

The repository contains multiple related analyses rather than one end-to-end automated pipeline. Run each stage only after its stated inputs have been supplied. The original code uses machine-specific local paths; these are not portable without editing.

## A. Full-domain data and nine-model forecasts

1. Acquire the original GLDAS-Noah V2.1 and GLDAS-DA V2.2 inputs from their providers. The base data-processing scripts use monthly GLDAS-Noah inputs and GLDAS-DA groundwater-storage information.
2. Edit raw-data paths, project roots and output directories in each relevant script.
3. Run `step1_extract_noah_data.m`, `step2_upscale_noah_data.m`, `step3_extract_DA_data.m`, `step4_upscale_DA_data.m`, and `step5_create_mask.m`, in that order. The final master database is expected at `<projectRoot>/s1_master_database/03_master_database/GLDAS_common_valid_mask_1deg_200401_202412.mat` for the additional experiments.
4. Run AR, ARIMA, SARIMA, SVM, RF, LSTM, SVMX, RFX, and LSTMX from the three model directories. Read the headers to set original/repeat modes, seed and parallel options; the 2022–2024 forecast tests use training/validation/test periods that advance without future predictors.
5. After saving all model prediction MAT files, run the scripts in `s2_analysis/` in stage order, checking each script's required paths and upstream MAT outputs. `s6_spatial_block_bootstrap.m` uses saved outputs and does not retrain models.

## B. Sensitivity and regional experiments

### `s6_SKLGB_sensitivity_2000_fixed_dualPC`: hyperparameter perturbations

- The provided `SKLGB_hyperparameter_sample2000.mat` fixes the 2,000-cell sample. PC1 takes alternating odd sample positions; PC2 takes even positions, 1,000 each. Do not resample.
- On **each** PC, set `projectRoot` near the start of `SKLGB_pilot10_run.m` to the full-domain original master database location.
- PC1: run `SKLGB_PC1_preflight` and then `SKLGB_PC1_train`.
- PC2: run `SKLGB_PC2_preflight` and then `SKLGB_PC2_train`.
- Each PC evaluates six one-at-a-time variants for each of RF, RFX, LSTM and LSTMX (24 perturbed configurations), across test years 2022–2024. The original all-model baseline is reused rather than globally retrained.
- After both PCs finish, transfer their `compact_sensitivity_results.mat` files to `merge_inputs/PC1/` and `merge_inputs/PC2/` respectively, then run `SKLGB_merge_2000_results` in MATLAB. The **two compact files already included here** allow this merge/summary stage to be rerun directly without repeating training. Generated `merged_results/` files have been omitted from the package.
- The `SKLGB_hyperparameter_sample400.mat` file is kept for an optional nested-sample consistency check during preflight. Do not change `block_num` after checkpoints exist.

### `s7_SKLGB_extend_2019_2021_singlePC`: extended early test years

- Edit `cfg.projectRoot` in `SKLGB_ext_config.m`.
- In MATLAB, `cd` to this folder, or add it to the MATLAB path. Run:

```matlab
SKLGB_ext_run('preflight');
SKLGB_ext_smoke('AR');
SKLGB_ext_smoke('ARIMA');
SKLGB_ext_smoke('SARIMA');
SKLGB_ext_run('train');
SKLGB_ext_run('analyze');
```

- This is a single-PC, fixed-2,000-cell run of all nine configurations for test years 2019, 2020, 2021. Main outputs appear in the `s7_extended_2019_2021_2000grids_singlePC` analysis directory under `cfg.projectRoot`. Source-specific AR/ARIMA/SARIMA core functions are included; old `Original_linear_sources/` reference copies were removed because execution does not use them.

### `s8_SKLGB_monthlymean_target_2000_singlePC`: monthly-mean target sensitivity

- Edit `cfg.projectRoot` and **`cfg.rawDAFolder`** in `SKLGB_mean_config.m`. This analysis needs the daily GLDAS-DA `GWS_tavg` files from 2003-12-01 to 2024-12-31, with 7,702 daily files under year directories, as well as the original master database and fixed sample.
- Run preflight and preprocessing, then smoke-test a few models before training:

```matlab
SKLGB_mean_run('preflight');
SKLGB_mean_run('preprocess');
SKLGB_mean_smoke('AR');
SKLGB_mean_smoke('RF');
SKLGB_mean_run('train');
SKLGB_mean_run('analyze');
```

- `preprocess` carries out native daily-to-monthly extraction, spherical-area-weighted 1-degree aggregation, and construction of a fixed-grid alternative-target master. It creates 253 monthly states and can consume substantial disk space (roughly 1–2 GB in the original experiment). All nine models use their original configurations, re-evaluated for 2022–2024. The original full-domain master and forecasts are not overwritten.

### `s9_Koppen_Climate_Analysis_Package`: Koppen-Geiger climate zones

- This stage performs regional summaries only; it does not train models.
- Supply `s3_uc_analysis1_results.mat` from the main `s2_analysis` workflow (not bundled), with its expected 2022–2024 arrays, and use the supplied `Koppen_1991_2020/koppen_geiger_1p0.tif`.
- Open `koppen_climate_attribution.m` and run. If the input MAT is not found locally, MATLAB prompts for it. Output CSV/optional XLSX/MAT/figures are written to `outputs/` alongside the script; this relative output path is the only portability adjustment applied to its original MATLAB source.

## Limitations of the release

- Input data and 17,284-cell original model outputs are not fully bundled. Reproducing all training requires substantial computing resources.
- Code was statically inspected and the package was validated as a ZIP, but native MATLAB execution was not performed in this environment. Read the individual source headers before running.
- Before distributing a newly generated `*.mat` file or model output, ensure it does not contain unintended local machine-specific metadata or restricted third-party data.
