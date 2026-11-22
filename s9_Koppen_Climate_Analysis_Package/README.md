# Regional attribution by Koppen-Geiger climate zone

This script computes climate-zone summaries for the original 2022–2024 nine-model evaluation without retraining. It uses five major Köppen-Geiger zones (A Tropical, B Arid, C Temperate, D Cold and E Polar), checks a fixed 17,284-cell global mask, and compares annual as well as three-year-mean attribution results.

## Included

- `koppen_climate_attribution.m` — main MATLAB script.
- `Koppen_1991_2020/koppen_geiger_1p0.tif` — 1-degree climate class map used in the analysis.

## Not included but required

- `s3_uc_analysis1_results.mat`, produced by the main `s2_analysis` pipeline. The script expects arrays for the original 2022–2024 evaluation, including `annualNRMSEStack`, `annualResult`, `commonTemporalMask`, coordinates, model names and test years. It will ask you to select this file if it is not found next to the script.

## Run

1. Put `s3_uc_analysis1_results.mat` next to this script or allow the MATLAB file-selection dialog to locate it.
2. Keep the supplied GeoTIFF in its `Koppen_1991_2020` subfolder.
3. Run `koppen_climate_attribution.m` in MATLAB.
4. Check the CSV files, optional XLSX file and figures created in this folder's `outputs/` subfolder. The output path was changed from the original computer-specific drive to this relative folder; the calculations are otherwise unchanged.

The script verifies the original global category counts, annual attribution consistency and gain-component reconstruction. Its zone-level summaries are conditional on the GLDAS valid-grid mask and use equal weights per retained cell, not an independently reweighted estimate of land area.

**Dataset reference:** Beck et al. (2023), *High-resolution (1 km) Koppen-Geiger maps for 1901-2099 based on constrained CMIP6 projections*, *Scientific Data*, 10:724. DOI: 10.1038/s41597-023-02549-6.
