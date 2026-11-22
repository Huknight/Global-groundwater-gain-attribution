%% s2_gain_attribution1.m

clc;
clearvars;

%
% ==================== 1. Code Description ====================
% Purpose:
%   Calculate and plot pixel-wise all-model-fusion results for:
%       1. Original run
%       2. UC1 run
%       3. UC2 run
%
% All-model fusion:
%   At each valid pixel, the model with the minimum NRMSE is selected from:
%       AR, ARIMA, SARIMA,
%       SVM, RF, LSTM,
%       SVMX, RFX, LSTMX.
%
% Figures:
%   Fig1 : Global best NRMSE of all-model fusion - Original run.
%   Fig2 : Global best NRMSE of all-model fusion - UC1 run.
%   Fig3 : Global best NRMSE of all-model fusion - UC2 run.
%
%   Fig4 : Global best-model class - Original run.
%   Fig5 : Global best-model class - UC1 run.
%   Fig6 : Global best-model class - UC2 run.
%
%   Fig7 : Best-model percentage pie chart - Original run.
%   Fig8 : Best-model percentage pie chart - UC1 run.
%   Fig9 : Best-model percentage pie chart - UC2 run.
%
%   Fig10: Pixel-wise standard deviation of the all-model-fusion best NRMSE
%          across the original, UC1, and UC2 runs.
%
%   Fig11: Pixel-wise stability of the best-model class across three runs.
%          Unchanged = the same best-model class in all three runs.
%          Changed   = at least one run has a different best-model class.
%
%   Fig12: Percentage pie chart for changed and unchanged best-model class.
%
% Outputs:
%   1. Twelve PNG figures.
%   2. MAT file containing all calculated maps and statistics.
%   3. CSV and Excel statistics tables.
%
% Data source:
%   All model results are read from the revised master database
%   (M_manuscript1_revised\s1_master_database), in which the 0.25-degree
%   to 1-degree upscaling has been corrected with spherical-area weights.
%
% Notes:
%   1. AR, ARIMA, and SARIMA are deterministic and use the same result
%      files in the original, UC1, and UC2 comparisons.
%   2. The six stochastic models use their corresponding original, UC1,
%      and UC2 result files.
%   3. A pixel is retained only when all nine models contain finite NRMSE
%      values for the corresponding run.
%   4. Fig10-Fig12 use only pixels valid in all three all-model-fusion runs.
%   5. If models have exactly identical NRMSE values, MATLAB min() selects
%      the model appearing first in the defined model order.
%
% ==================== 2. Centralized Parameter Settings ====================

% ------------------------------ Paths ------------------------------
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';

dataRoot = fullfile(projectRoot, 's1_master_database');

maskFile = fullfile(dataRoot, ...
    '03_master_database', ...
    'GLDAS_common_valid_mask_1deg_200401_202412.mat');

linearDir = fullfile(dataRoot, '04_linear_models');
nonlinearDir = fullfile(dataRoot, '05_nonlinear_models');
nonlinearXDir = fullfile(dataRoot, '06_nonlinear_x_models');

outDir = fullfile(projectRoot, ...
    's3_darft', ...
    's1_analysis', ...
    's2_gain_attribution1_results');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ------------------ Linear-model result files ------------------
% These deterministic models are shared by all three runs.
% Note: the revised SARIMA result file has no "_d0" suffix.

arFile = fullfile(linearDir, ...
    'AR_univariate_global_200401_202412.mat');

arimaFile = fullfile(linearDir, ...
    'ARIMA_univariate_global_200401_202412_d0.mat');

sarimaFile = fullfile(linearDir, ...
    'SARIMA_univariate_global_200401_202412.mat');

% ------------------ Original-run result files ------------------
svmFile = fullfile(nonlinearDir, ...
    'SVM_univariate_global_200401_202412.mat');

rfFile = fullfile(nonlinearDir, ...
    'RF_univariate_global_200401_202412.mat');

lstmFile = fullfile(nonlinearDir, ...
    'LSTM_univariate_global_200401_202412.mat');

svmxFile = fullfile(nonlinearXDir, ...
    'SVMX_multivariate_global_200401_202412.mat');

rfxFile = fullfile(nonlinearXDir, ...
    'RFX_multivariate_global_200401_202412.mat');

lstmxFile = fullfile(nonlinearXDir, ...
    'LSTMX_multivariate_global_200401_202412.mat');

% ------------------------ UC1 result files ------------------------
svmUC1File = fullfile(nonlinearDir, ...
    'SVM_univariate_global_200401_202412_uc1.mat');

rfUC1File = fullfile(nonlinearDir, ...
    'RF_univariate_global_200401_202412_uc1.mat');

lstmUC1File = fullfile(nonlinearDir, ...
    'LSTM_univariate_global_200401_202412_uc1.mat');

svmxUC1File = fullfile(nonlinearXDir, ...
    'SVMX_multivariate_global_200401_202412_uc1.mat');

rfxUC1File = fullfile(nonlinearXDir, ...
    'RFX_multivariate_global_200401_202412_uc1.mat');

lstmxUC1File = fullfile(nonlinearXDir, ...
    'LSTMX_multivariate_global_200401_202412_uc1.mat');

% ------------------------ UC2 result files ------------------------
svmUC2File = fullfile(nonlinearDir, ...
    'SVM_univariate_global_200401_202412_uc2.mat');

rfUC2File = fullfile(nonlinearDir, ...
    'RF_univariate_global_200401_202412_uc2.mat');

lstmUC2File = fullfile(nonlinearDir, ...
    'LSTM_univariate_global_200401_202412_uc2.mat');

svmxUC2File = fullfile(nonlinearXDir, ...
    'SVMX_multivariate_global_200401_202412_uc2.mat');

rfxUC2File = fullfile(nonlinearXDir, ...
    'RFX_multivariate_global_200401_202412_uc2.mat');

lstmxUC2File = fullfile(nonlinearXDir, ...
    'LSTMX_multivariate_global_200401_202412_uc2.mat');

% ------------------ Model order and variable names ------------------
modelNames = { ...
    'AR', ...
    'ARIMA', ...
    'SARIMA', ...
    'SVM', ...
    'RF', ...
    'LSTM', ...
    'SVMX', ...
    'RFX', ...
    'LSTMX'};

nrmseVarNames = { ...
    'AR_NRMSE_mean', ...
    'ARIMA_NRMSE_mean', ...
    'SARIMA_NRMSE_mean', ...
    'SVM_NRMSE_mean', ...
    'RF_NRMSE_mean', ...
    'LSTM_NRMSE_mean', ...
    'SVMX_NRMSE_mean', ...
    'RFX_NRMSE_mean', ...
    'LSTMX_NRMSE_mean'};

nModel = numel(modelNames);

% ----------------------- Result-file groups -----------------------

% ---------- Original run ----------
originalFiles = { ...
    arFile, arimaFile, sarimaFile, ...
    svmFile, rfFile, lstmFile, ...
    svmxFile, rfxFile, lstmxFile};

% ---------- UC1 run ----------
uc1Files = { ...
    arFile, arimaFile, sarimaFile, ...
    svmUC1File, rfUC1File, lstmUC1File, ...
    svmxUC1File, rfxUC1File, lstmxUC1File};

% ---------- UC2 run ----------
uc2Files = { ...
    arFile, arimaFile, sarimaFile, ...
    svmUC2File, rfUC2File, lstmUC2File, ...
    svmxUC2File, rfxUC2File, lstmxUC2File};

% ----------------------- Shared model colors -----------------------

% ---------- Linear models: red shades ----------
cAR = [1.00, 0.75, 0.75];
cARIMA = [1.00, 0.30, 0.30];
cSARIMA = [0.60, 0.00, 0.00];

% ---------- Nonlinear models: blue shades ----------
cSVM = [0.70, 0.85, 1.00];
cRF = [0.20, 0.55, 1.00];
cLSTM = [0.00, 0.20, 0.70];

% ---------- Nonlinear-X models: green shades ----------
cSVMX = [0.75, 1.00, 0.75];
cRFX = [0.30, 0.75, 0.30];
cLSTMX = [0.00, 0.45, 0.00];

all9Cmap = [ ...
    cAR;
    cARIMA;
    cSARIMA;
    cSVM;
    cRF;
    cLSTM;
    cSVMX;
    cRFX;
    cLSTMX];

all9Labels = modelNames;

% ------------------- Stability-class colors -------------------

% Class 1: unchanged best-model type across three runs
% Class 2: changed best-model type across three runs

stabilityLabels = { ...
    'Unchanged', ...
    'Changed'};

stabilityCmap = [ ...
    0.00, 0.45, 0.90;   % blue: unchanged
    0.85, 0.10, 0.10];  % red: changed

% ---------------------- Unified figure settings ----------------------
settings.mapPosition = [100, 100, 1000, 700];
settings.piePosition = [200, 200, 100, 100];

settings.bestNRMSEClim = [0, 0.4];
settings.fusionSTDClim = [0, 0.02];

settings.bestNRMSECmap = jet(256);
settings.fusionSTDCmap = jet(256);

settings.fontName = 'Times New Roman';

settings.titleFontSize = 16;
settings.axisLabelFontSize = 13;
settings.tickFontSize = 13;
settings.colorbarFontSize = 13;
settings.colorbarLabelFontSize = 13;

settings.pieFontName = 'Arial';
settings.pieFontSize = 6;

settings.axisLineWidth = 0.8;
settings.pieEdgeLineWidth = 0.5;

settings.outputResolution = 300;

%
% ==================== 3. Main Body ====================

% ------------------- Load the common valid mask -------------------
D_mask = load(maskFile, 'CommonValid');
validMask = logical(D_mask.CommonValid);
clear D_mask

% ------------------- Load latitude and longitude -------------------
% Coordinates are read from the AR result file.
D_coord = load(arFile, 'lat', 'lon');
lat = double(D_coord.lat(:));
lon = double(D_coord.lon(:));
clear D_coord

[Lon2D, Lat2D] = meshgrid(lon, lat);

% ---------------- Calculate original-run all-model fusion ----------------
fprintf('\nCalculating original-run all-model fusion...\n');

[allFusionOriginal, bestModelOriginal, completeMaskOriginal] = ...
    calculateAllModelFusion(originalFiles, nrmseVarNames, validMask);

% --------------------- Calculate UC1 all-model fusion ---------------------
fprintf('Calculating UC1 all-model fusion...\n');

[allFusionUC1, bestModelUC1, completeMaskUC1] = ...
    calculateAllModelFusion(uc1Files, nrmseVarNames, validMask);

% --------------------- Calculate UC2 all-model fusion ---------------------
fprintf('Calculating UC2 all-model fusion...\n');

[allFusionUC2, bestModelUC2, completeMaskUC2] = ...
    calculateAllModelFusion(uc2Files, nrmseVarNames, validMask);

% ------------------- Fig1: Original-run best NRMSE -------------------
fig1File = fullfile(outDir, 'Fig1_AllModelFusion_BestNRMSE_Original.png');

plotContinuousMap( ...
    Lon2D, Lat2D, allFusionOriginal, ...
    'Global best NRMSE of all-model fusion - Original run', ...
    'Best NRMSE', ...
    settings.bestNRMSEClim, settings.bestNRMSECmap, ...
    fig1File, settings);

% ----------------------- Fig2: UC1 best NRMSE -----------------------
fig2File = fullfile(outDir, 'Fig2_AllModelFusion_BestNRMSE_UC1.png');

plotContinuousMap( ...
    Lon2D, Lat2D, allFusionUC1, ...
    'Global best NRMSE of all-model fusion - UC1 run', ...
    'Best NRMSE', ...
    settings.bestNRMSEClim, settings.bestNRMSECmap, ...
    fig2File, settings);

% ----------------------- Fig3: UC2 best NRMSE -----------------------
fig3File = fullfile(outDir, 'Fig3_AllModelFusion_BestNRMSE_UC2.png');

plotContinuousMap( ...
    Lon2D, Lat2D, allFusionUC2, ...
    'Global best NRMSE of all-model fusion - UC2 run', ...
    'Best NRMSE', ...
    settings.bestNRMSEClim, settings.bestNRMSECmap, ...
    fig3File, settings);

% ------------- Fig4: Original-run best-model class -------------
fig4File = fullfile(outDir, 'Fig4_AllModelFusion_BestModelClass_Original.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, bestModelOriginal, ...
    all9Cmap, all9Labels, ...
    'Global best model class of all-model fusion - Original run', ...
    'Best model', [0.5, 9.5], ...
    fig4File, settings);

% ----------------- Fig5: UC1 best-model class -----------------
fig5File = fullfile(outDir, 'Fig5_AllModelFusion_BestModelClass_UC1.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, bestModelUC1, ...
    all9Cmap, all9Labels, ...
    'Global best model class of all-model fusion - UC1 run', ...
    'Best model', [0.5, 9.5], ...
    fig5File, settings);

% ----------------- Fig6: UC2 best-model class -----------------
fig6File = fullfile(outDir, 'Fig6_AllModelFusion_BestModelClass_UC2.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, bestModelUC2, ...
    all9Cmap, all9Labels, ...
    'Global best model class of all-model fusion - UC2 run', ...
    'Best model', [0.5, 9.5], ...
    fig6File, settings);

% ------------------- Calculate best-model counts -------------------
bestModelCountsOriginal = countClasses( ...
    bestModelOriginal, completeMaskOriginal, nModel);

bestModelCountsUC1 = countClasses( ...
    bestModelUC1, completeMaskUC1, nModel);

bestModelCountsUC2 = countClasses( ...
    bestModelUC2, completeMaskUC2, nModel);

% ---------------- Calculate best-model percentages ----------------
bestModelPercentOriginal = calculatePercentages(bestModelCountsOriginal);
bestModelPercentUC1 = calculatePercentages(bestModelCountsUC1);
bestModelPercentUC2 = calculatePercentages(bestModelCountsUC2);

% ------------- Fig7: Original-run best-model pie chart -------------
fig7File = fullfile(outDir, 'Fig7_AllModelFusion_BestModelPercentage_Original.png');

plotPiePercent(bestModelCountsOriginal, all9Cmap, fig7File, settings);

% ----------------- Fig8: UC1 best-model pie chart -----------------
fig8File = fullfile(outDir, 'Fig8_AllModelFusion_BestModelPercentage_UC1.png');

plotPiePercent(bestModelCountsUC1, all9Cmap, fig8File, settings);

% ----------------- Fig9: UC2 best-model pie chart -----------------
fig9File = fullfile(outDir, 'Fig9_AllModelFusion_BestModelPercentage_UC2.png');

plotPiePercent(bestModelCountsUC2, all9Cmap, fig9File, settings);

% ------------------- Common three-run valid mask -------------------

% Fig10-Fig12 use only pixels valid in all three all-model-fusion runs.

commonThreeRunMask = ...
    completeMaskOriginal & ...
    completeMaskUC1 & ...
    completeMaskUC2 & ...
    isfinite(allFusionOriginal) & ...
    isfinite(allFusionUC1) & ...
    isfinite(allFusionUC2) & ...
    isfinite(bestModelOriginal) & ...
    isfinite(bestModelUC1) & ...
    isfinite(bestModelUC2);

% ------------------- Fig10: Best-NRMSE uncertainty -------------------

% Stack the best NRMSE maps from the three independent runs.
allFusionStack = cat(3, allFusionOriginal, allFusionUC1, allFusionUC2);

% Calculate pixel-wise sample standard deviation using N-1.
allFusionSTD = std(allFusionStack, 0, 3);

allFusionSTD(~commonThreeRunMask) = NaN;

fig10File = fullfile(outDir, 'Fig10_AllModelFusion_BestNRMSE_STD.png');

plotContinuousMap( ...
    Lon2D, Lat2D, allFusionSTD, ...
    'Global best-NRMSE uncertainty of all-model fusion', ...
    'Best NRMSE standard deviation', ...
    settings.fusionSTDClim, settings.fusionSTDCmap, ...
    fig10File, settings);

% ------------------- Fig11: Best-model stability class -------------------

% A pixel is unchanged only when the exact best-model class is identical
% in the original, UC1, and UC2 runs.

unchangedBestModelMask = ...
    commonThreeRunMask & ...
    (bestModelOriginal == bestModelUC1) & ...
    (bestModelOriginal == bestModelUC2);

changedBestModelMask = ...
    commonThreeRunMask & ...
    ~unchangedBestModelMask;

% Class definition:
%   1 = unchanged
%   2 = changed
bestModelStabilityMap = nan(size(validMask));

bestModelStabilityMap(unchangedBestModelMask) = 1;
bestModelStabilityMap(changedBestModelMask) = 2;

fig11File = fullfile(outDir, 'Fig11_AllModelFusion_BestModel_Stability.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, bestModelStabilityMap, ...
    stabilityCmap, stabilityLabels, ...
    'Three-run stability of all-model-fusion best-model class', ...
    'Best-model stability', [0.5, 2.5], ...
    fig11File, settings);

% --------------- Fig12: Best-model stability pie chart ---------------

stabilityCounts = [ ...
    nnz(unchangedBestModelMask), ...
    nnz(changedBestModelMask)];

stabilityPercentages = calculatePercentages(stabilityCounts);

fig12File = fullfile(outDir, 'Fig12_AllModelFusion_BestModel_Stability_Percentage.png');

plotPiePercent(stabilityCounts, stabilityCmap, fig12File, settings);

% ------------------- Best-model statistics table -------------------
bestModelTable = table( ...
    modelNames(:), ...
    bestModelCountsOriginal(:), ...
    bestModelPercentOriginal(:), ...
    bestModelCountsUC1(:), ...
    bestModelPercentUC1(:), ...
    bestModelCountsUC2(:), ...
    bestModelPercentUC2(:), ...
    'VariableNames', { ...
    'Model', ...
    'Original_Count', ...
    'Original_Percent', ...
    'UC1_Count', ...
    'UC1_Percent', ...
    'UC2_Count', ...
    'UC2_Percent'});

% ------------------- Stability statistics table -------------------
stabilityTable = table( ...
    stabilityLabels(:), ...
    stabilityCounts(:), ...
    stabilityPercentages(:), ...
    'VariableNames', { ...
    'BestModelStability', ...
    'PixelCount', ...
    'PixelPercent'});

% ------------------- Best-NRMSE uncertainty statistics -------------------
validFusionSTD = allFusionSTD(commonThreeRunMask);

meanFusionSTD = mean(validFusionSTD);
medianFusionSTD = median(validFusionSTD);
percentile10FusionSTD = prctile(validFusionSTD, 10);
percentile90FusionSTD = prctile(validFusionSTD, 90);

fusionSTDTable = table( ...
    meanFusionSTD, ...
    medianFusionSTD, ...
    percentile10FusionSTD, ...
    percentile90FusionSTD, ...
    nnz(commonThreeRunMask), ...
    'VariableNames', { ...
    'Mean_BestNRMSE_STD', ...
    'Median_BestNRMSE_STD', ...
    'P10_BestNRMSE_STD', ...
    'P90_BestNRMSE_STD', ...
    'ValidPixelCount'});

% ------------------- Display command-window statistics -------------------
fprintf('\n');
fprintf('==========================================================================\n');
fprintf('                 ALL-MODEL-FUSION BEST-MODEL STATISTICS\n');
fprintf('==========================================================================\n');

fprintf('%-8s %10s %10s %10s %10s %10s %10s\n', ...
    'Model', ...
    'Ori Count', ...
    'Ori (%)', ...
    'UC1 Count', ...
    'UC1 (%)', ...
    'UC2 Count', ...
    'UC2 (%)');

fprintf('--------------------------------------------------------------------------\n');

for iModel = 1:nModel
    fprintf('%-8s %10d %10.1f %10d %10.1f %10d %10.1f\n', ...
        modelNames{iModel}, ...
        bestModelCountsOriginal(iModel), ...
        bestModelPercentOriginal(iModel), ...
        bestModelCountsUC1(iModel), ...
        bestModelPercentUC1(iModel), ...
        bestModelCountsUC2(iModel), ...
        bestModelPercentUC2(iModel));
end

fprintf('==========================================================================\n');

fprintf('\nValid all-model pixels:\n');
fprintf('Original run : %d\n', nnz(completeMaskOriginal));
fprintf('UC1 run      : %d\n', nnz(completeMaskUC1));
fprintf('UC2 run      : %d\n', nnz(completeMaskUC2));
fprintf('All 3 runs   : %d\n', nnz(commonThreeRunMask));

fprintf('\nGlobal mean all-model-fusion NRMSE:\n');

fprintf('Original run : %.6f\n', ...
    mean(allFusionOriginal(completeMaskOriginal), 'omitnan'));

fprintf('UC1 run      : %.6f\n', ...
    mean(allFusionUC1(completeMaskUC1), 'omitnan'));

fprintf('UC2 run      : %.6f\n', ...
    mean(allFusionUC2(completeMaskUC2), 'omitnan'));

fprintf('\nBest-NRMSE uncertainty across three runs:\n');
fprintf('Mean STD     : %.6f\n', meanFusionSTD);
fprintf('Median STD   : %.6f\n', medianFusionSTD);
fprintf('P10 STD      : %.6f\n', percentile10FusionSTD);
fprintf('P90 STD      : %.6f\n', percentile90FusionSTD);

fprintf('\nBest-model class stability across three runs:\n');

fprintf('Unchanged    : %d pixels (%.2f%%)\n', ...
    stabilityCounts(1), stabilityPercentages(1));

fprintf('Changed      : %d pixels (%.2f%%)\n', ...
    stabilityCounts(2), stabilityPercentages(2));

disp(' ');
disp(bestModelTable);

disp(' ');
disp(stabilityTable);

disp(' ');
disp(fusionSTDTable);

% ------------------- Save CSV tables -------------------
bestModelCsvFile = fullfile(outDir, 's2_gain_attribution1_best_model_statistics.csv');
writetable(bestModelTable, bestModelCsvFile);

stabilityCsvFile = fullfile(outDir, 's2_gain_attribution1_best_model_stability.csv');
writetable(stabilityTable, stabilityCsvFile);

fusionSTDCsvFile = fullfile(outDir, 's2_gain_attribution1_best_NRMSE_STD_statistics.csv');
writetable(fusionSTDTable, fusionSTDCsvFile);

% ------------------- Save Excel tables -------------------
excelFile = fullfile(outDir, 's2_gain_attribution1_statistics.xlsx');

writetable(bestModelTable, excelFile, 'Sheet', 'BestModel');
writetable(stabilityTable, excelFile, 'Sheet', 'ModelStability');
writetable(fusionSTDTable, excelFile, 'Sheet', 'BestNRMSE_STD');

% ------------------- Save calculation results -------------------
matFile = fullfile(outDir, 's2_gain_attribution1_results.mat');

save(matFile, ...
    'modelNames', ...
    'nrmseVarNames', ...
    'all9Cmap', ...
    'all9Labels', ...
    'stabilityCmap', ...
    'stabilityLabels', ...
    'allFusionOriginal', ...
    'allFusionUC1', ...
    'allFusionUC2', ...
    'allFusionSTD', ...
    'allFusionStack', ...
    'bestModelOriginal', ...
    'bestModelUC1', ...
    'bestModelUC2', ...
    'bestModelStabilityMap', ...
    'completeMaskOriginal', ...
    'completeMaskUC1', ...
    'completeMaskUC2', ...
    'commonThreeRunMask', ...
    'unchangedBestModelMask', ...
    'changedBestModelMask', ...
    'bestModelCountsOriginal', ...
    'bestModelCountsUC1', ...
    'bestModelCountsUC2', ...
    'bestModelPercentOriginal', ...
    'bestModelPercentUC1', ...
    'bestModelPercentUC2', ...
    'stabilityCounts', ...
    'stabilityPercentages', ...
    'meanFusionSTD', ...
    'medianFusionSTD', ...
    'percentile10FusionSTD', ...
    'percentile90FusionSTD', ...
    'bestModelTable', ...
    'stabilityTable', ...
    'fusionSTDTable', ...
    'validMask', ...
    'lat', ...
    'lon', ...
    'settings');

fprintf('\nCompleted: all-model-fusion performance and stability analysis.\n');
fprintf('Output directory:\n%s\n', outDir);

%
% ==================== 4. Additional Functions ====================

function [bestNRMSE, bestModelIdx, completeModelMask] = ...
    calculateAllModelFusion(modelFiles, nrmseVarNames, validMask)
% At each pixel, select the model with the minimum NRMSE from a given group
% of model result files. A pixel is complete only when all models contain
% finite NRMSE values.

    nModel = numel(modelFiles);

    nRow = size(validMask, 1);
    nCol = size(validMask, 2);

    modelStack = nan(nRow, nCol, nModel);

    %% ---------- Load all model NRMSE maps ----------

    for iModel = 1:nModel
        thisFile = modelFiles{iModel};
        thisVariable = nrmseVarNames{iModel};

        D_model = load(thisFile, thisVariable);
        thisNRMSE = double(D_model.(thisVariable));

        modelStack(:, :, iModel) = thisNRMSE;

        clear D_model thisNRMSE
    end

    %% ---------- Identify pixels valid for all models ----------

    completeModelMask = validMask;

    for iModel = 1:nModel
        completeModelMask = ...
            completeModelMask & ...
            isfinite(modelStack(:, :, iModel));
    end

    %% ---------- Select the minimum NRMSE and model index ----------

    [bestNRMSE, bestModelIdx] = min(modelStack, [], 3);

    bestNRMSE(~completeModelMask) = NaN;
    bestModelIdx(~completeModelMask) = NaN;
end

function plotContinuousMap( ...
    Lon2D, Lat2D, mapData, ...
    figureTitle, colorbarLabel, colorLimits, cmapUse, ...
    figureFile, settings)
% Plot a continuous map with pcolor and export it to a PNG file.

    %% ---------- Create figure ----------

    fig = figure('Color', 'w', 'Position', settings.mapPosition);

    ax = axes('Parent', fig);

    %% ---------- Plot continuous map ----------

    hMap = pcolor(ax, Lon2D, Lat2D, mapData);

    set(hMap, 'EdgeColor', 'none');

    %% ---------- Axis settings ----------

    set(ax, ...
        'YDir', 'normal', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.tickFontSize, ...
        'LineWidth', settings.axisLineWidth);

    axis(ax, 'equal');
    axis(ax, 'tight');

    xlabel(ax, 'Longitude', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    ylabel(ax, 'Latitude', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    title(ax, figureTitle, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.titleFontSize, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none');

    %% ---------- Colormap and limits ----------

    colormap(ax, cmapUse);
    clim(ax, colorLimits);

    %% ---------- Colorbar ----------

    cb = colorbar(ax);
    cb.FontName = settings.fontName;
    cb.FontSize = settings.colorbarFontSize;

    ylabel(cb, colorbarLabel, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.colorbarLabelFontSize);

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', settings.outputResolution);
end

function plotDiscreteMap( ...
    Lon2D, Lat2D, mapData, ...
    cmapUse, classLabels, ...
    figureTitle, colorbarLabel, colorLimits, ...
    figureFile, settings)
% Plot a discrete class map with pcolor and export it to a PNG file.

    %% ---------- Create figure ----------

    fig = figure('Color', 'w', 'Position', settings.mapPosition);

    ax = axes('Parent', fig);

    %% ---------- Plot discrete map ----------

    hMap = pcolor(ax, Lon2D, Lat2D, mapData);

    set(hMap, 'EdgeColor', 'none');

    %% ---------- Axis settings ----------

    set(ax, ...
        'YDir', 'normal', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.tickFontSize, ...
        'LineWidth', settings.axisLineWidth);

    axis(ax, 'equal');
    axis(ax, 'tight');

    xlabel(ax, 'Longitude', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    ylabel(ax, 'Latitude', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    title(ax, figureTitle, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.titleFontSize, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none');

    %% ---------- Discrete colors ----------

    colormap(ax, cmapUse);
    clim(ax, colorLimits);

    %% ---------- Colorbar ----------

    cb = colorbar(ax);

    cb.Ticks = 1:numel(classLabels);
    cb.TickLabels = classLabels;
    cb.FontName = settings.fontName;
    cb.FontSize = settings.colorbarFontSize;

    ylabel(cb, colorbarLabel, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.colorbarLabelFontSize);

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', settings.outputResolution);
end

function counts = countClasses(classMap, validClassMask, nClass)
% Count the number of valid pixels in each discrete class.

    counts = zeros(1, nClass);

    for iClass = 1:nClass
        counts(iClass) = nnz(validClassMask & classMap == iClass);
    end
end

function percentages = calculatePercentages(counts)
% Convert class counts to percentages of the total count.

    counts = double(counts(:)');

    totalCount = sum(counts);

    if totalCount <= 0
        percentages = nan(size(counts));
    else
        percentages = counts ./ totalCount .* 100;
    end
end

function plotPiePercent(counts, cmapUse, figureFile, settings)
% Plot a percentage pie chart of class counts and export it to a PNG file.

    counts = double(counts(:)');

    totalCount = sum(counts);

    if totalCount <= 0
        warning(['The pie chart was not created because all counts are zero.']);
        return;
    end

    %% ---------- Percentage labels ----------

    percentages = counts ./ totalCount .* 100;

    labels = cell(size(counts));

    for iClass = 1:numel(counts)
        labels{iClass} = sprintf('%.1f%%', percentages(iClass));
    end

    %% ---------- Create figure ----------

    fig = figure('Color', 'w', 'Position', settings.piePosition);

    ax = axes('Parent', fig);

    positiveIdx = find(counts > 0);

    if isempty(positiveIdx)
        close(fig);
        return;
    end

    %% ---------- Plot pie chart ----------

    h = pie(ax, counts, labels);

    axis(ax, 'equal');
    axis(ax, 'off');

    patchHandles = h(1:2:end);
    textHandles = h(2:2:end);

    %% ---------- Set slice colors ----------

    for iPatch = 1:numel(patchHandles)
        thisClassIdx = positiveIdx(iPatch);

        set(patchHandles(iPatch), ...
            'FaceColor', cmapUse(thisClassIdx, :), ...
            'EdgeColor', 'w', ...
            'LineWidth', settings.pieEdgeLineWidth);
    end

    %% ---------- Set text properties ----------

    for iText = 1:numel(textHandles)
        set(textHandles(iText), ...
            'FontSize', settings.pieFontSize, ...
            'FontName', settings.pieFontName, ...
            'Color', 'k');
    end

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', settings.outputResolution);
end
