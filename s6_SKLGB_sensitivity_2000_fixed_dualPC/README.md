# Fixed 2,000-grid hyperparameter sensitivity (two PCs)

This is the experiment reported for key RF/RFX and LSTM/LSTMX hyperparameters, using a **fixed 2,000-cell spatially stratified sample** and an immutable two-computer partition. PC1 evaluates sample positions 1,3,...,1999; PC2 positions 2,4,...,2000. Each PC fits its 1,000 cells over the 2022, 2023, 2024 test years.

## Required inputs

- `SKLGB_hyperparameter_sample2000.mat`: fixed sample and reference per-cell metadata (included).
- `SKLGB_hyperparameter_sample400.mat`: retained for optional nested-sample audit (included).
- `SKLGB_global_original_reference.mat`: original full-domain summary for anchored comparison (included).
- Full preprocessed GLDAS master at `<projectRoot>/s1_master_database/03_master_database/GLDAS_common_valid_mask_1deg_200401_202412.mat` (not included).
- MATLAB toolboxes for bagged regression trees, deep learning and parallel computing when enabled.

## Starting a new two-PC run

1. Copy the folder and identical sample files to both computers.
2. Edit **only the machine-specific `projectRoot`** and, if necessary, worker count near the top of `SKLGB_pilot10_run.m`; keep `nPilot=2000`, `machineCount=2`, `mode='perturb_only'`.
3. On PC1 run `SKLGB_PC1_preflight`, verify `PREFLIGHT PASSED`, then `SKLGB_PC1_train`.
4. On PC2 run `SKLGB_PC2_preflight`, verify `PREFLIGHT PASSED`, then `SKLGB_PC2_train`.
5. Do not alter the sample, `block_num`, settings or partial checkpoints when resuming the same experiment.

For each of RF and RFX, six variants modify tree count or minimum leaf size; for LSTM and LSTMX, six variants modify hidden units or initial learning rate. Each variant repeats the original validation-based lag search. The full model-fitting process is expensive.

## Reproduce merge and summary with bundled compact outputs

The subdirectories `merge_inputs/PC1/` and `merge_inputs/PC2/` contain the **two already completed compact result MAT files**, one per machine. From MATLAB with this directory on the path, run:

```matlab
SKLGB_merge_2000_results
```

The merger checks each PC's disjoint cell identity, fixed sample membership, configurations, and completed test-year records. It generates a new `merged_results/` directory with model ranking, gain category, attribution components, runtime reports and consistency records. `SKLGB_analyze_merged_2000.m` is invoked by the merge program. Previously generated `merged_results/` files were removed from this source package because they can be regenerated.

The 2,000-grid results are stratified-sample estimates and should not be described as new full-domain retraining of all 17,284 cells.
