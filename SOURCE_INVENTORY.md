# Source inventory and release notes

- MATLAB source files: **79** (29 from the base archive; 50 retained from the four added folders).
- Supporting data files: **8** (fixed sample/reference/compact merge inputs and one climate GeoTIFF).
- New English-only operating guides: root README/RUNNING plus one README in each added experiment folder.
- No MATLAB model-fitting, statistical, or attribution algorithms were changed.
- The only executable source edits were to `s9_Koppen_Climate_Analysis_Package/koppen_climate_attribution.m`: match the real file name and GeoTIFF folder in comments, and use a relative `outputs/` path instead of a Windows D: output location.
- Omitted: generated `merged_results/` files, run logs, 0-byte PC readmes, Chinese/old README .txt files, ancillary coordinate/allocation CSVs, and archival reference copies under `Original_linear_sources/` or `Original_processing_references/`.
- Although the two compact s6 result MATs are included for re-aggregation, the complete global forecast outputs and raw GLDAS products are not bundled.

## Retained MATLAB source files

| File | SHA-256 (first 16 characters) |
|---|---|
| `s1_scripts/s1_data_prepare/step1_extract_noah_data.m` | `3011a81f96b3d575` |
| `s1_scripts/s1_data_prepare/step2_upscale_noah_data.m` | `a86101e3f174c8c5` |
| `s1_scripts/s1_data_prepare/step3_extract_DA_data.m` | `cad51d674093d96f` |
| `s1_scripts/s1_data_prepare/step4_upscale_DA_data.m` | `eefc3347d7712401` |
| `s1_scripts/s1_data_prepare/step5_create_mask.m` | `d34bc2dfcfcedccf` |
| `s1_scripts/s2_linear_models/ARIMA_global_forecast.m` | `bf6d6dcab7579b09` |
| `s1_scripts/s2_linear_models/AR_global_forecast.m` | `9e2e7aa03d31321a` |
| `s1_scripts/s2_linear_models/SARIMA_global_forecast.m` | `f86d3da4cb1441d1` |
| `s1_scripts/s3_nonlinear_models/LSTM_global_forecast.m` | `8997b7dc369a7173` |
| `s1_scripts/s3_nonlinear_models/RF_global_forecast.m` | `012db87a5e472bda` |
| `s1_scripts/s3_nonlinear_models/SVM_gs_global_forecast.m` | `7fbe86b664e5c666` |
| `s1_scripts/s4_nonlinear_x_models/LSTMX_global_forecast.m` | `2e9a58fefcec8f68` |
| `s1_scripts/s4_nonlinear_x_models/RFX_global_forecast.m` | `c048f3f2badedd41` |
| `s1_scripts/s4_nonlinear_x_models/SVMX_global_forecast.m` | `3e20deb8dd1b359d` |
| `s2_analysis/s1_overall_performance1.m` | `28c327f44032a7f0` |
| `s2_analysis/s1_overall_performance2.m` | `06e0bb5f5fbfa2fb` |
| `s2_analysis/s1_overall_performance3.m` | `fc2b73306362b5f6` |
| `s2_analysis/s2_gain_attribution1.m` | `1f4ee06081c15d50` |
| `s2_analysis/s2_gain_attribution2.m` | `115ccdc50a7957c4` |
| `s2_analysis/s2_gain_attribution3.m` | `c811ff817917ec7a` |
| `s2_analysis/s3_uc_analysis1.m` | `4c789fc9894d88d4` |
| `s2_analysis/s3_uc_analysis1_1.m` | `3836996ca84fefd1` |
| `s2_analysis/s3_uc_analysis2.m` | `ab7e98b641069b30` |
| `s2_analysis/s4_NRMSE_selection.m` | `73a709c11b18096d` |
| `s2_analysis/s4_model_selection.m` | `5b33a31701e50c68` |
| `s2_analysis/s5_NRMSE_sensitivity1.m` | `078907d74b210678` |
| `s2_analysis/s5_NRMSE_sensitivity2.m` | `90d219b86f02f4fb` |
| `s2_analysis/s5_NRMSE_sensitivity3.m` | `e01314a37eb5d6f9` |
| `s2_analysis/s6_spatial_block_bootstrap.m` | `d57976929a1d3b7e` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_PC1_preflight.m` | `ba99bba28490fa8a` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_PC1_train.m` | `de7bc6b646639392` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_PC2_preflight.m` | `7a8a6acdce9101f0` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_PC2_train.m` | `1606d4da8c194943` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_analyze_merged_2000.m` | `211b76538645ff5d` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_export_compact_results.m` | `f0dba32beafcb3f0` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_merge_2000_results.m` | `366e04253704a4dc` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot10_prepare_master.m` | `632539e69ff32fb9` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot10_run.m` | `4c3a3a0fd80ea184` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot10_summary.m` | `7aaf0150aa690aa3` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot_LSTM.m` | `4fd8aa529f86bc9a` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot_LSTMX.m` | `e54aed3a8417f43a` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot_RF.m` | `5915d1eaf01e6166` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_pilot_RFX.m` | `d22379ed99691a74` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_plot_sampling_map.m` | `6cad8f50b78f9ac2` |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_validate_selection.m` | `04fac4b2b43bd666` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_ARIMA_core.m` | `d451b010f2659a7f` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_AR_core.m` | `9fbd12d370328ee2` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_LSTM.m` | `e456eaa43538949e` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_LSTMX.m` | `e16740a2be30df57` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_RF.m` | `64db1004cd3f7bed` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_RFX.m` | `57ae6bdfec7a98b5` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_SARIMA_core.m` | `db1ba8729145514b` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_SVM.m` | `e1c79b52b39620cf` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_SVMX.m` | `509e49f6cb3b2639` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_analyze.m` | `414efc3dc857098e` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_config.m` | `8c78f32c8d69abfa` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_linear.m` | `c97ed920f8f03c2b` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_period_bootstrap.m` | `3ecbc2b584d3ac9f` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_prepare_master.m` | `e4d220f4b94caa01` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_run.m` | `a2597c01d783a34e` |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_ext_smoke.m` | `54500267ac523a14` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_ext_ARIMA_core.m` | `d451b010f2659a7f` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_ext_AR_core.m` | `9fbd12d370328ee2` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_ext_SARIMA_core.m` | `db1ba8729145514b` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_LSTM.m` | `7463a59b35992caa` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_LSTMX.m` | `cb97cb8e4c611630` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_RF.m` | `0c72718c332d1aa0` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_RFX.m` | `94a1c882870c5c03` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_SVM.m` | `2984a35326da7d1a` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_SVMX.m` | `33f70bc597d912d1` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_analyze.m` | `aa9bc121f3f42bb4` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_config.m` | `2696ebbef6f9b5f9` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_extract.m` | `69341964923e2b3b` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_linear.m` | `b182327ef8c9954b` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_prepare.m` | `7fb2fa31bdbfcefc` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_run.m` | `0d138941e8ae175f` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_smoke.m` | `c9970cd77bf02ac6` |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_mean_upscale.m` | `d1eb0255528fc121` |
| `s9_Koppen_Climate_Analysis_Package/koppen_climate_attribution.m` | `b32e6bc7a9a42ac7` |

## Retained supporting files

| File | Size (KiB) |
|---|---:|
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_global_original_reference.mat` | 2.4 |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_hyperparameter_sample2000.mat` | 440.2 |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/SKLGB_hyperparameter_sample400.mat` | 90.6 |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/merge_inputs/PC1/compact_sensitivity_results.mat` | 1532.1 |
| `s6_SKLGB_sensitivity_2000_fixed_dualPC/merge_inputs/PC2/compact_sensitivity_results.mat` | 1532.4 |
| `s7_SKLGB_extend_2019_2021_singlePC/SKLGB_hyperparameter_sample2000.mat` | 440.2 |
| `s8_SKLGB_monthlymean_target_2000_singlePC/SKLGB_hyperparameter_sample2000.mat` | 440.2 |
| `s9_Koppen_Climate_Analysis_Package/Koppen_1991_2020/koppen_geiger_1p0.tif` | 9.3 |

## Portability edit audit

`s9_Koppen_Climate_Analysis_Package/koppen_climate_attribution.m`:

- Original SHA-256: `25ae7f7555ef210867e693b44d29dbd921061dfc4c50d0b88424da95a18354a2`
- Packaged SHA-256: `b32e6bc7a9a42ac7f654868dd823d85dd3143a1afb288747dce0b25dc32b02ab`
