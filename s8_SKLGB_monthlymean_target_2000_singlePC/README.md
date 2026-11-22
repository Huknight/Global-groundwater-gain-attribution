# Alternative target: monthly-mean storage change (one PC)

This sensitivity experiment replaces the monthly-end-derived groundwater-storage-change target with a target defined by differences between consecutive **1-degree monthly-mean GLDAS-DA groundwater storage states**, while using exactly the same fixed 2,000 valid cells and original nine model configurations for 2022–2024.

## Inputs

- The fixed sample `SKLGB_hyperparameter_sample2000.mat` (included).
- The original full-domain GLDAS master MAT (not included).
- GLDAS-DA daily NetCDF files containing `GWS_tavg` for 2003-12-01 to 2024-12-31 (**7,702 files**, as expected by the original workflow); the files must be arranged under year folders and are not included.
- MATLAB toolboxes consistent with the nine original models.

## Workflow

1. Edit the project root, raw daily GLDAS directory and optional parallel settings in `SKLGB_mean_config.m`. The original example `cfg.rawDAFolder` contains a computer-specific path and must be replaced.
2. From this folder in MATLAB, run:

```matlab
SKLGB_mean_run('preflight');
SKLGB_mean_run('preprocess');
```

Preflight checks the original reference sample, native inputs and the daily input files. Preprocess builds monthly mean states, 1-degree spatial aggregates and a new 2,000-cell alternative-target master. The month-by-month extraction is resumable and saves 253 monthly state files; expect substantial disk requirements (approximately 1–2 GB in the original run).

3. Before training, run smoke tests as needed, for example:

```matlab
SKLGB_mean_smoke('AR');
SKLGB_mean_smoke('RF');
```

4. Train all nine configurations and compare with the original month-end target:

```matlab
SKLGB_mean_run('train');
SKLGB_mean_run('analyze');
```

Output is under `<projectRoot>/s3_darft/s1_analysis/s8_deltaGWS_monthmean_2000grids_singlePC/`. The original target master and existing global runs are not overwritten. The `Original_processing_references/` folder of historical script copies was removed; the active extraction/upscaling/model routines are preserved under this folder.
