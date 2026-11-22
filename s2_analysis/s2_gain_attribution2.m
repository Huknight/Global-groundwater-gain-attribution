%% s2_gain_attribution2.m

clc;
clearvars;

%
% ==================== 1. Code Description ====================
% Purpose:
%   Calculate and plot global gain-attribution categories for:
%       1. Original run
%       2. UC1 run
%       3. UC2 run
%
% Input:
%   1. s2_gain_attribution1_results.mat
%      Used for coordinates, masks, and previous best-model results.
%
%   2. Nine model NRMSE result files
%      Required to reconstruct group-level winner and second-best logic.
%
% Attribution classes:
%   1. No gain
%   2. Linear gain
%   3. Nonlinear gain
%   4. Hydroclimatic gain
%   5. Mixed gain
%
% Group definitions:
%   Linear:
%       AR, ARIMA, SARIMA
%
%   Nonlinear:
%       SVM, RF, LSTM
%
%   NonlinearX:
%       SVMX, RFX, LSTMX
%
% Attribution logic:
%   No gain:
%       Best group = Linear and best linear model = AR.
%
%   Linear gain:
%       Best group = Linear and best linear model = ARIMA or SARIMA.
%
%   Nonlinear gain:
%       Best group = Nonlinear.
%
%   Hydroclimatic gain:
%       Best group = NonlinearX;
%       second-best group = Nonlinear;
%       best NonlinearX model is paired with the best Nonlinear model.
%
%   Mixed gain:
%       Best group = NonlinearX, but the clean hydroclimatic conditions
%       are not satisfied.
%
% Figures:
%   Fig1: Global gain-category distribution - Original run.
%   Fig2: Global gain-category distribution - UC1 run.
%   Fig3: Global gain-category distribution - UC2 run.
%
%   Fig4: Gain-category percentage pie chart - Original run.
%   Fig5: Gain-category percentage pie chart - UC1 run.
%   Fig6: Gain-category percentage pie chart - UC2 run.
%
%   Fig7: Three-run gain-category stability map.
%         Blue = unchanged across all three runs.
%         Red  = changed in at least one run.
%
%   Fig8: Gain-category stability percentage pie chart.
%
% Outputs:
%   1. Eight PNG figures.
%   2. CSV and Excel summary tables.
%   3. MAT file containing all calculated maps and masks.
%
% Data source:
%   All model results are read from the revised master database
%   (M_manuscript1_revised\s1_master_database), in which the 0.25-degree
%   to 1-degree upscaling has been corrected with spherical-area weights.
%   The previous all-model-fusion results are read from the revised
%   s2_gain_attribution1_results.mat produced by s2_gain_attribution1.m
%   on the same revised database.
%
% Notes:
%   1. Strict stability is used:
%      Original category = UC1 category = UC2 category.
%
%   2. Mixed-L and Mixed-N are combined into one Mixed-gain class.
%
%   3. The script verifies that the reconstructed all-model winner agrees
%      with the winner saved by s2_gain_attribution1.m.
%
%   4. If two models have exactly identical NRMSE values, MATLAB min()
%      selects the model appearing first in the defined model order.
%
% ==================== 2. Centralized Parameter Settings ====================

% ------------------------------ Paths ------------------------------
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';

dataRoot = fullfile(projectRoot, 's1_master_database');

linearDir = fullfile(dataRoot, '04_linear_models');
nonlinearDir = fullfile(dataRoot, '05_nonlinear_models');
nonlinearXDir = fullfile(dataRoot, '06_nonlinear_x_models');

% Previous all-model-fusion results from the revised database
previousResultDir = fullfile(projectRoot, ...
    's3_darft', ...
    's1_analysis', ...
    's2_gain_attribution1_results');

previousResultFile = fullfile(previousResultDir, ...
    's2_gain_attribution1_results.mat');

outDir = fullfile(projectRoot, ...
    's3_darft', ...
    's1_analysis', ...
    's2_gain_attribution2_results');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ------------------ Linear-model result files ------------------
% These deterministic models are shared across all three runs.
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

% ------------------- Gain-category definitions -------------------
gainClassNames = { ...
    'No gain', ...
    'Linear gain', ...
    'Nonlinear gain', ...
    'Hydroclimatic gain', ...
    'Mixed gain'};

gainClassShortNames = { ...
    'No gain', ...
    'Linear', ...
    'Nonlinear', ...
    'Hydroclimatic', ...
    'Mixed'};

nGainClass = numel(gainClassNames);

% ------------------- Gain-category colors -------------------

% Class 1: No gain
% Class 2: Linear gain
% Class 3: Nonlinear gain
% Class 4: Hydroclimatic gain
% Class 5: Mixed gain

gainClassCmap = [ ...
    0.75, 0.75, 0.75;   % No gain: light gray
    0.85, 0.10, 0.10;   % Linear gain: red
    0.00, 0.35, 0.90;   % Nonlinear gain: blue
    0.00, 0.60, 0.20;   % Hydroclimatic gain: green
    1.00, 0.80, 0.00];  % Mixed gain: yellow

% ------------------- Stability-class definitions -------------------
stabilityLabels = { ...
    'Unchanged', ...
    'Changed'};

stabilityCmap = [ ...
    0.00, 0.45, 0.90;   % Unchanged: blue
    0.85, 0.10, 0.10];  % Changed: red

% ---------------------- Unified figure settings ----------------------
settings.mapPosition = [100, 100, 1000, 700];
settings.piePosition = [200, 200, 100, 100];

settings.fontName = 'Times New Roman';

settings.titleFontSize = 15;
settings.axisLabelFontSize = 13;
settings.tickFontSize = 11;
settings.colorbarFontSize = 11;
settings.colorbarLabelFontSize = 12;

settings.pieFontName = 'Arial';
settings.pieFontSize = 6;

settings.axisLineWidth = 0.8;
settings.pieEdgeLineWidth = 0.5;

settings.outputResolution = 300;

%
% ==================== 3. Main Body ====================

% ------------------- Load previous attribution-1 results -------------------

% Variables loaded from the revised attribution-1 results:
%   validMask, lat, lon,
%   bestModelOriginal, bestModelUC1, bestModelUC2,
%   completeMaskOriginal, completeMaskUC1, completeMaskUC2

requiredPreviousVariables = { ...
    'validMask', ...
    'lat', ...
    'lon', ...
    'bestModelOriginal', ...
    'bestModelUC1', ...
    'bestModelUC2', ...
    'completeMaskOriginal', ...
    'completeMaskUC1', ...
    'completeMaskUC2'};

S_previous = load(previousResultFile, requiredPreviousVariables{:});

validMask = logical(S_previous.validMask);

lat = double(S_previous.lat(:));
lon = double(S_previous.lon(:));

previousBestModelOriginal = double(S_previous.bestModelOriginal);
previousBestModelUC1 = double(S_previous.bestModelUC1);
previousBestModelUC2 = double(S_previous.bestModelUC2);

completeMaskOriginal = logical(S_previous.completeMaskOriginal);
completeMaskUC1 = logical(S_previous.completeMaskUC1);
completeMaskUC2 = logical(S_previous.completeMaskUC2);

clear S_previous

[Lon2D, Lat2D] = meshgrid(lon, lat);

% ------------------- Original-run attribution -------------------
fprintf('\nCalculating original-run gain attribution...\n');

R_original = calculateGainAttribution( ...
    originalFiles, nrmseVarNames, validMask, completeMaskOriginal);

% ----------------------- UC1 attribution -----------------------
fprintf('Calculating UC1 gain attribution...\n');

R_uc1 = calculateGainAttribution( ...
    uc1Files, nrmseVarNames, validMask, completeMaskUC1);

% ----------------------- UC2 attribution -----------------------
fprintf('Calculating UC2 gain attribution...\n');

R_uc2 = calculateGainAttribution( ...
    uc2Files, nrmseVarNames, validMask, completeMaskUC2);

% ------------------- Check reconstructed winners -------------------
winnerMismatchOriginal = nnz( ...
    completeMaskOriginal & ...
    R_original.bestAll9Idx ~= previousBestModelOriginal);

winnerMismatchUC1 = nnz( ...
    completeMaskUC1 & ...
    R_uc1.bestAll9Idx ~= previousBestModelUC1);

winnerMismatchUC2 = nnz( ...
    completeMaskUC2 & ...
    R_uc2.bestAll9Idx ~= previousBestModelUC2);

if winnerMismatchOriginal > 0
    warning(['Original-run reconstructed winner differs from the previous ', ...
        'result at %d pixels.'], winnerMismatchOriginal);
end

if winnerMismatchUC1 > 0
    warning(['UC1 reconstructed winner differs from the previous result ', ...
        'at %d pixels.'], winnerMismatchUC1);
end

if winnerMismatchUC2 > 0
    warning(['UC2 reconstructed winner differs from the previous result ', ...
        'at %d pixels.'], winnerMismatchUC2);
end

% ------------------- Fig1: Original-run gain categories -------------------
fig1File = fullfile(outDir, 'Fig1_Global_Gain_Category_Original.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, R_original.gainClassMap, ...
    gainClassCmap, gainClassShortNames, ...
    'Global gain-category distribution - Original run', ...
    'Gain category', [0.5, 5.5], ...
    fig1File, settings);

% ----------------------- Fig2: UC1 gain categories -----------------------
fig2File = fullfile(outDir, 'Fig2_Global_Gain_Category_UC1.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, R_uc1.gainClassMap, ...
    gainClassCmap, gainClassShortNames, ...
    'Global gain-category distribution - UC1 run', ...
    'Gain category', [0.5, 5.5], ...
    fig2File, settings);

% ----------------------- Fig3: UC2 gain categories -----------------------
fig3File = fullfile(outDir, 'Fig3_Global_Gain_Category_UC2.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, R_uc2.gainClassMap, ...
    gainClassCmap, gainClassShortNames, ...
    'Global gain-category distribution - UC2 run', ...
    'Gain category', [0.5, 5.5], ...
    fig3File, settings);

% ------------------- Gain-category counts -------------------
gainCountsOriginal = countClasses( ...
    R_original.gainClassMap, completeMaskOriginal, nGainClass);

gainCountsUC1 = countClasses( ...
    R_uc1.gainClassMap, completeMaskUC1, nGainClass);

gainCountsUC2 = countClasses( ...
    R_uc2.gainClassMap, completeMaskUC2, nGainClass);

% ------------------- Gain-category percentages -------------------
gainPercentOriginal = calculatePercentages(gainCountsOriginal);
gainPercentUC1 = calculatePercentages(gainCountsUC1);
gainPercentUC2 = calculatePercentages(gainCountsUC2);

% ------------------- Fig4: Original-run gain pie chart -------------------
fig4File = fullfile(outDir, 'Fig4_Global_Gain_Category_Percentage_Original.png');

plotPiePercent(gainCountsOriginal, gainClassCmap, fig4File, settings);

% ----------------------- Fig5: UC1 gain pie chart -----------------------
fig5File = fullfile(outDir, 'Fig5_Global_Gain_Category_Percentage_UC1.png');

plotPiePercent(gainCountsUC1, gainClassCmap, fig5File, settings);

% ----------------------- Fig6: UC2 gain pie chart -----------------------
fig6File = fullfile(outDir, 'Fig6_Global_Gain_Category_Percentage_UC2.png');

plotPiePercent(gainCountsUC2, gainClassCmap, fig6File, settings);

% ------------------- Three-run common attribution mask -------------------
commonThreeRunMask = ...
    completeMaskOriginal & ...
    completeMaskUC1 & ...
    completeMaskUC2 & ...
    isfinite(R_original.gainClassMap) & ...
    isfinite(R_uc1.gainClassMap) & ...
    isfinite(R_uc2.gainClassMap);

% ------------------- Fig7: Gain-category stability -------------------

% A category is unchanged only when all three runs have exactly the same
% attribution-class number.

unchangedGainClassMask = ...
    commonThreeRunMask & ...
    R_original.gainClassMap == R_uc1.gainClassMap & ...
    R_original.gainClassMap == R_uc2.gainClassMap;

changedGainClassMask = ...
    commonThreeRunMask & ...
    ~unchangedGainClassMask;

% Class definition:
%   1 = unchanged
%   2 = changed
gainClassStabilityMap = nan(size(validMask));

gainClassStabilityMap(unchangedGainClassMask) = 1;
gainClassStabilityMap(changedGainClassMask) = 2;

fig7File = fullfile(outDir, 'Fig7_ThreeRun_Gain_Category_Stability.png');

plotDiscreteMap( ...
    Lon2D, Lat2D, gainClassStabilityMap, ...
    stabilityCmap, stabilityLabels, ...
    'Three-run stability of all-model-fusion gain category', ...
    'Gain-category stability', [0.5, 2.5], ...
    fig7File, settings);

% ------------------- Fig8: Gain-category stability pie chart -------------------
stabilityCounts = [ ...
    nnz(unchangedGainClassMask), ...
    nnz(changedGainClassMask)];

stabilityPercentages = calculatePercentages(stabilityCounts);

fig8File = fullfile(outDir, 'Fig8_ThreeRun_Gain_Category_Stability_Percentage.png');

plotPiePercent(stabilityCounts, stabilityCmap, fig8File, settings);

% ------------------- Gain-category summary table -------------------
gainCategoryTable = table( ...
    gainClassNames(:), ...
    gainCountsOriginal(:), ...
    gainPercentOriginal(:), ...
    gainCountsUC1(:), ...
    gainPercentUC1(:), ...
    gainCountsUC2(:), ...
    gainPercentUC2(:), ...
    'VariableNames', { ...
    'GainCategory', ...
    'Original_Count', ...
    'Original_Percent', ...
    'UC1_Count', ...
    'UC1_Percent', ...
    'UC2_Count', ...
    'UC2_Percent'});

% ------------------- Stability summary table -------------------
stabilityTable = table( ...
    stabilityLabels(:), ...
    stabilityCounts(:), ...
    stabilityPercentages(:), ...
    'VariableNames', { ...
    'GainCategoryStability', ...
    'PixelCount', ...
    'PixelPercent'});

% ------------------- Attribution diagnostics table -------------------
runNames = { ...
    'Original'; ...
    'UC1'; ...
    'UC2'};

validPixelCount = [ ...
    nnz(completeMaskOriginal); ...
    nnz(completeMaskUC1); ...
    nnz(completeMaskUC2)];

assignedPixelCount = [ ...
    R_original.nAssigned; ...
    R_uc1.nAssigned; ...
    R_uc2.nAssigned];

unassignedPixelCount = [ ...
    R_original.nUnassigned; ...
    R_uc1.nUnassigned; ...
    R_uc2.nUnassigned];

multiAssignedPixelCount = [ ...
    R_original.nMultiAssigned; ...
    R_uc1.nMultiAssigned; ...
    R_uc2.nMultiAssigned];

winnerMismatchCount = [ ...
    winnerMismatchOriginal; ...
    winnerMismatchUC1; ...
    winnerMismatchUC2];

diagnosticsTable = table( ...
    runNames, ...
    validPixelCount, ...
    assignedPixelCount, ...
    unassignedPixelCount, ...
    multiAssignedPixelCount, ...
    winnerMismatchCount, ...
    'VariableNames', { ...
    'Run', ...
    'ValidPixels', ...
    'AssignedPixels', ...
    'UnassignedPixels', ...
    'MultiAssignedPixels', ...
    'WinnerMismatchPixels'});

% ------------------- Display command-window results -------------------
fprintf('\n');
fprintf('============================================================================\n');
fprintf('                    GAIN-ATTRIBUTION CATEGORY STATISTICS\n');
fprintf('============================================================================\n');

fprintf('%-22s %10s %10s %10s %10s %10s %10s\n', ...
    'Gain category', ...
    'Ori Count', ...
    'Ori (%)', ...
    'UC1 Count', ...
    'UC1 (%)', ...
    'UC2 Count', ...
    'UC2 (%)');

fprintf('----------------------------------------------------------------------------\n');

for iClass = 1:nGainClass
    fprintf('%-22s %10d %10.2f %10d %10.2f %10d %10.2f\n', ...
        gainClassNames{iClass}, ...
        gainCountsOriginal(iClass), ...
        gainPercentOriginal(iClass), ...
        gainCountsUC1(iClass), ...
        gainPercentUC1(iClass), ...
        gainCountsUC2(iClass), ...
        gainPercentUC2(iClass));
end

fprintf('============================================================================\n');

fprintf('\nThree-run gain-category stability:\n');

fprintf('Unchanged : %d pixels (%.2f%%)\n', ...
    stabilityCounts(1), stabilityPercentages(1));

fprintf('Changed   : %d pixels (%.2f%%)\n', ...
    stabilityCounts(2), stabilityPercentages(2));

fprintf('\nWinner reconstruction checks:\n');

fprintf('Original winner mismatches : %d\n', winnerMismatchOriginal);
fprintf('UC1 winner mismatches      : %d\n', winnerMismatchUC1);
fprintf('UC2 winner mismatches      : %d\n', winnerMismatchUC2);

fprintf('\nAttribution assignment checks:\n');

fprintf('Original: assigned = %d, unassigned = %d, multi-assigned = %d\n', ...
    R_original.nAssigned, R_original.nUnassigned, R_original.nMultiAssigned);

fprintf('UC1     : assigned = %d, unassigned = %d, multi-assigned = %d\n', ...
    R_uc1.nAssigned, R_uc1.nUnassigned, R_uc1.nMultiAssigned);

fprintf('UC2     : assigned = %d, unassigned = %d, multi-assigned = %d\n', ...
    R_uc2.nAssigned, R_uc2.nUnassigned, R_uc2.nMultiAssigned);

disp(' ');
disp(gainCategoryTable);

disp(' ');
disp(stabilityTable);

disp(' ');
disp(diagnosticsTable);

% ------------------- Save CSV tables -------------------
gainCategoryCsvFile = fullfile(outDir, 's2_gain_attribution2_gain_category_statistics.csv');
writetable(gainCategoryTable, gainCategoryCsvFile);

stabilityCsvFile = fullfile(outDir, 's2_gain_attribution2_gain_category_stability.csv');
writetable(stabilityTable, stabilityCsvFile);

diagnosticsCsvFile = fullfile(outDir, 's2_gain_attribution2_diagnostics.csv');
writetable(diagnosticsTable, diagnosticsCsvFile);

% ------------------- Save Excel workbook -------------------
excelFile = fullfile(outDir, 's2_gain_attribution2_statistics.xlsx');

writetable(gainCategoryTable, excelFile, 'Sheet', 'GainCategories');
writetable(stabilityTable, excelFile, 'Sheet', 'Stability');
writetable(diagnosticsTable, excelFile, 'Sheet', 'Diagnostics');

% ------------------- Save calculation results -------------------
matFile = fullfile(outDir, 's2_gain_attribution2_results.mat');

save(matFile, ...
    'modelNames', ...
    'nrmseVarNames', ...
    'gainClassNames', ...
    'gainClassShortNames', ...
    'gainClassCmap', ...
    'stabilityLabels', ...
    'stabilityCmap', ...
    'R_original', ...
    'R_uc1', ...
    'R_uc2', ...
    'gainCountsOriginal', ...
    'gainCountsUC1', ...
    'gainCountsUC2', ...
    'gainPercentOriginal', ...
    'gainPercentUC1', ...
    'gainPercentUC2', ...
    'commonThreeRunMask', ...
    'unchangedGainClassMask', ...
    'changedGainClassMask', ...
    'gainClassStabilityMap', ...
    'stabilityCounts', ...
    'stabilityPercentages', ...
    'gainCategoryTable', ...
    'stabilityTable', ...
    'diagnosticsTable', ...
    'validMask', ...
    'completeMaskOriginal', ...
    'completeMaskUC1', ...
    'completeMaskUC2', ...
    'lat', ...
    'lon', ...
    'settings', ...
    '-v7.3');

fprintf('\nCompleted: three-run gain-attribution category analysis.\n');
fprintf('Output directory:\n%s\n', outDir);

%
% ==================== 4. Additional Functions ====================

function R = calculateGainAttribution( ...
    modelFiles, nrmseVarNames, validMask, completeModelMask)
% Reconstruct the group-level winners and assign each pixel to one of the
% five gain-attribution classes.

    nModel = numel(modelFiles);

    nRow = size(validMask, 1);
    nCol = size(validMask, 2);

    modelStack = nan(nRow, nCol, nModel);

    %% ---------- Load nine NRMSE maps ----------

    for iModel = 1:nModel
        thisFile = modelFiles{iModel};
        thisVariable = nrmseVarNames{iModel};

        D_model = load(thisFile, thisVariable);
        thisNRMSE = double(D_model.(thisVariable));

        modelStack(:, :, iModel) = thisNRMSE;

        clear D_model thisNRMSE
    end

    %% ---------- Confirm complete-model mask ----------

    calculatedCompleteMask = validMask;

    for iModel = 1:nModel
        calculatedCompleteMask = ...
            calculatedCompleteMask & ...
            isfinite(modelStack(:, :, iModel));
    end

    analysisMask = ...
        validMask & ...
        completeModelMask & ...
        calculatedCompleteMask;

    %% ---------- Construct group stacks ----------

    linearStack = modelStack(:, :, 1:3);
    nonlinearStack = modelStack(:, :, 4:6);
    nonlinearXStack = modelStack(:, :, 7:9);

    %% ---------- Calculate group-level fusion ----------

    [linearFusion, bestLinearLocalIdx] = min(linearStack, [], 3);
    [nonlinearFusion, bestNonlinearLocalIdx] = min(nonlinearStack, [], 3);
    [nonlinearXFusion, bestNonlinearXLocalIdx] = min(nonlinearXStack, [], 3);

    linearFusion(~analysisMask) = NaN;
    nonlinearFusion(~analysisMask) = NaN;
    nonlinearXFusion(~analysisMask) = NaN;

    bestLinearLocalIdx(~analysisMask) = NaN;
    bestNonlinearLocalIdx(~analysisMask) = NaN;
    bestNonlinearXLocalIdx(~analysisMask) = NaN;

    %% ---------- Convert to global model indices ----------

    bestLinearGlobalIdx = bestLinearLocalIdx;
    bestNonlinearGlobalIdx = bestNonlinearLocalIdx + 3;
    bestNonlinearXGlobalIdx = bestNonlinearXLocalIdx + 6;

    bestLinearGlobalIdx(~analysisMask) = NaN;
    bestNonlinearGlobalIdx(~analysisMask) = NaN;
    bestNonlinearXGlobalIdx(~analysisMask) = NaN;

    %% ---------- Rank the three model groups ----------

    groupFusionStack = cat(3, linearFusion, nonlinearFusion, nonlinearXFusion);

    [sortedGroupNRMSE, sortedGroupIdx] = sort(groupFusionStack, 3, 'ascend');

    bestGroupIdx = sortedGroupIdx(:, :, 1);
    secondGroupIdx = sortedGroupIdx(:, :, 2);

    bestGroupNRMSE = sortedGroupNRMSE(:, :, 1);

    bestGroupIdx(~analysisMask) = NaN;
    secondGroupIdx(~analysisMask) = NaN;
    bestGroupNRMSE(~analysisMask) = NaN;

    %% ---------- Reconstruct all-model winner ----------

    bestAll9Idx = nan(size(validMask));

    maskBestLinear = analysisMask & bestGroupIdx == 1;
    maskBestNonlinear = analysisMask & bestGroupIdx == 2;
    maskBestNonlinearX = analysisMask & bestGroupIdx == 3;

    bestAll9Idx(maskBestLinear) = bestLinearGlobalIdx(maskBestLinear);
    bestAll9Idx(maskBestNonlinear) = bestNonlinearGlobalIdx(maskBestNonlinear);
    bestAll9Idx(maskBestNonlinearX) = bestNonlinearXGlobalIdx(maskBestNonlinearX);

    allFusion = bestGroupNRMSE;
    allFusion(~analysisMask) = NaN;

    %% ---------- Paired hydroclimatic logic ----------

    % Local pairing:
    %   1: SVM  <-> SVMX
    %   2: RF   <-> RFX
    %   3: LSTM <-> LSTMX

    maskHydroClean = ...
        analysisMask & ...
        bestGroupIdx == 3 & ...
        secondGroupIdx == 2 & ...
        bestNonlinearXLocalIdx == bestNonlinearLocalIdx;

    %% ---------- Gain-category masks ----------

    % Class 1: No gain
    maskNoGain = ...
        analysisMask & ...
        bestGroupIdx == 1 & ...
        bestLinearLocalIdx == 1;

    % Class 2: Linear gain
    maskLinearGain = ...
        analysisMask & ...
        bestGroupIdx == 1 & ...
        (bestLinearLocalIdx == 2 | bestLinearLocalIdx == 3);

    % Class 3: Nonlinear gain
    maskNonlinearGain = ...
        analysisMask & ...
        bestGroupIdx == 2;

    % Class 4: Hydroclimatic gain
    % Already defined above as maskHydroClean.

    % Class 5: Mixed gain
    % Combines:
    %   1. NonlinearX best + Linear second-best
    %   2. NonlinearX best + Nonlinear second-best + unpaired algorithms
    maskMixedGain = ...
        analysisMask & ...
        bestGroupIdx == 3 & ...
        ~maskHydroClean;

    %% ---------- Create five-class attribution map ----------

    gainClassMap = nan(size(validMask));

    gainClassMap(maskNoGain) = 1;
    gainClassMap(maskLinearGain) = 2;
    gainClassMap(maskNonlinearGain) = 3;
    gainClassMap(maskHydroClean) = 4;
    gainClassMap(maskMixedGain) = 5;

    %% ---------- Assignment checks ----------

    classSum = ...
        double(maskNoGain) + ...
        double(maskLinearGain) + ...
        double(maskNonlinearGain) + ...
        double(maskHydroClean) + ...
        double(maskMixedGain);

    nAssigned = nnz(analysisMask & classSum == 1);
    nUnassigned = nnz(analysisMask & classSum == 0);
    nMultiAssigned = nnz(analysisMask & classSum > 1);

    %% ---------- Store results ----------

    R.analysisMask = analysisMask;

    R.linearFusion = linearFusion;
    R.nonlinearFusion = nonlinearFusion;
    R.nonlinearXFusion = nonlinearXFusion;

    R.bestLinearLocalIdx = bestLinearLocalIdx;
    R.bestNonlinearLocalIdx = bestNonlinearLocalIdx;
    R.bestNonlinearXLocalIdx = bestNonlinearXLocalIdx;

    R.bestGroupIdx = bestGroupIdx;
    R.secondGroupIdx = secondGroupIdx;

    R.bestAll9Idx = bestAll9Idx;
    R.allFusion = allFusion;

    R.maskNoGain = maskNoGain;
    R.maskLinearGain = maskLinearGain;
    R.maskNonlinearGain = maskNonlinearGain;
    R.maskHydroClean = maskHydroClean;
    R.maskMixedGain = maskMixedGain;

    R.gainClassMap = gainClassMap;

    R.nAssigned = nAssigned;
    R.nUnassigned = nUnassigned;
    R.nMultiAssigned = nMultiAssigned;
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

    %% ---------- Colormap and limits ----------

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

    %% ---------- Set percentage-label properties ----------

    for iText = 1:numel(textHandles)
        set(textHandles(iText), ...
            'FontSize', settings.pieFontSize, ...
            'FontName', settings.pieFontName, ...
            'Color', 'k');
    end

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', settings.outputResolution);
end
