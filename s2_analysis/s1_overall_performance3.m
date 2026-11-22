%% s1_overall_performance3.m

clc;
clearvars;

%
% ==================== 1. Code Description ====================
% Purpose:
%   Quantify pixel-wise run-to-run variability in NRMSE for six stochastic
%   groundwater forecasting models using three independent model runs.
%
% Uncertainty metric:
%   Pixel-wise sample standard deviation of NRMSE across:
%       1. Original run
%       2. UC1 run
%       3. UC2 run
%
% Figures:
%   Fig1-Fig6: Global NRMSE-STD of SVM, RF, LSTM, SVMX, RFX, and LSTMX
%   models, respectively.
%
% Statistical outputs:
%   For each model, calculate the spatial:
%       1. Mean NRMSE standard deviation
%       2. Median NRMSE standard deviation
%       3. 10th percentile
%       4. 90th percentile
%
% Outputs:
%   1. Six global NRMSE-STD maps
%   2. One CSV statistics table
%   3. One MAT result file
%
% Data source:
%   All model results are read from the revised master database
%   (M_manuscript1_revised\s1_master_database), in which the 0.25-degree
%   to 1-degree upscaling has been corrected with spherical-area weights.
%
% Notes:
%   1. Only pixels with valid NRMSE values in all three runs are used.
%   2. A shared fixed color range [0, 0.04] is applied to all six maps so
%      that the maps are directly comparable.
%   3. All six maps use the jet colormap.
%   4. This analysis represents stochastic run-to-run variability rather
%      than complete predictive uncertainty.
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
    's1_overall_performance3_results');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ------------------- Coordinate reference file -------------------
% Latitude and longitude are loaded from the AR result file.
coordFile = fullfile(linearDir, ...
    'AR_univariate_global_200401_202412.mat');

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

% ------------------ Model names and variable names ------------------
modelNames = { ...
    'SVM', ...
    'RF', ...
    'LSTM', ...
    'SVMX', ...
    'RFX', ...
    'LSTMX'};

nrmseVarNames = { ...
    'SVM_NRMSE_mean', ...
    'RF_NRMSE_mean', ...
    'LSTM_NRMSE_mean', ...
    'SVMX_NRMSE_mean', ...
    'RFX_NRMSE_mean', ...
    'LSTMX_NRMSE_mean'};

% Each row contains:
%   column 1: original-run file
%   column 2: UC1 file
%   column 3: UC2 file
modelFiles = { ...
    svmFile,   svmUC1File,   svmUC2File;
    rfFile,    rfUC1File,    rfUC2File;
    lstmFile,  lstmUC1File,  lstmUC2File;
    svmxFile,  svmxUC1File,  svmxUC2File;
    rfxFile,   rfxUC1File,   rfxUC2File;
    lstmxFile, lstmxUC1File, lstmxUC2File};

nModel = numel(modelNames);
nRun = 3;

% ------------------------- Figure settings -------------------------
% Figure position and size
figPosition = [100, 100, 1000, 700];

% Shared colormap
uncertaintyColormap = jet(256);

% Font settings
fontName = 'Times New Roman';

titleFontSize = 16;
axisLabelFontSize = 13;
tickFontSize = 13;
colorbarFontSize = 13;
colorbarLabelFontSize = 13;

% Output resolution
outputResolution = 300;

%
% ==================== 3. Main Body ====================

% ------------------- Load the common valid mask -------------------
D_mask = load(maskFile, 'CommonValid');
validMask = logical(D_mask.CommonValid);
clear D_mask

% ------------------- Load latitude and longitude -------------------
% Coordinates are read from the AR result file.
D_coord = load(coordFile, 'lat', 'lon');
lat = double(D_coord.lat(:));
lon = double(D_coord.lon(:));
clear D_coord

% Construct two-dimensional longitude and latitude grids
[Lon2D, Lat2D] = meshgrid(lon, lat);

% ------------------- Initialize result arrays -------------------

% Pixel-wise NRMSE standard-deviation maps:
%   latitude x longitude x model
NRMSE_STD_maps = nan(size(validMask, 1), size(validMask, 2), nModel);

% Number of pixels with valid values in all three runs
validPixelCounts = zeros(nModel, 1);

% Spatial summary statistics
meanSTD = nan(nModel, 1);
medianSTD = nan(nModel, 1);
percentile10STD = nan(nModel, 1);
percentile90STD = nan(nModel, 1);

% ------------------- Calculate pixel-wise NRMSE uncertainty -------------------
for iModel = 1:nModel
    modelName = modelNames{iModel};
    nrmseVariable = nrmseVarNames{iModel};

    fprintf('\nProcessing %s model...\n', modelName);

    % Store the three run-specific NRMSE maps
    runStack = nan(size(validMask, 1), size(validMask, 2), nRun);

    for iRun = 1:nRun
        thisFile = modelFiles{iModel, iRun};

        D_model = load(thisFile, nrmseVariable);
        thisNRMSE = double(D_model.(nrmseVariable));

        runStack(:, :, iRun) = thisNRMSE;

        clear D_model thisNRMSE
    end

    %% ---------- Identify pixels valid in all three runs ----------

    completeRunMask = validMask & ...
        isfinite(runStack(:, :, 1)) & ...
        isfinite(runStack(:, :, 2)) & ...
        isfinite(runStack(:, :, 3));

    validPixelCounts(iModel) = nnz(completeRunMask);

    %% ---------- Calculate sample standard deviation ----------

    % The second input is 0, so MATLAB uses N-1 as the denominator.
    % The third input specifies that the calculation is performed across
    % the three independent runs.
    stdMap = std(runStack, 0, 3);

    % Retain only pixels valid in all three runs
    stdMap(~completeRunMask) = NaN;

    NRMSE_STD_maps(:, :, iModel) = stdMap;

    %% ---------- Calculate spatial summary statistics ----------

    validSTDValues = stdMap(completeRunMask);

    meanSTD(iModel) = mean(validSTDValues);
    medianSTD(iModel) = median(validSTDValues);
    percentile10STD(iModel) = prctile(validSTDValues, 10);
    percentile90STD(iModel) = prctile(validSTDValues, 90);

    clear runStack completeRunMask stdMap validSTDValues
end

% ------------------- Create uncertainty statistics table -------------------
uncertaintyTable = table( ...
    modelNames(:), ...
    meanSTD, ...
    medianSTD, ...
    percentile10STD, ...
    percentile90STD, ...
    'VariableNames', { ...
    'Model', ...
    'Mean_NRMSE_STD', ...
    'Median_NRMSE_STD', ...
    'P10_NRMSE_STD', ...
    'P90_NRMSE_STD'});

% ------------------- Display statistics in command window -------------------
fprintf('\n');
fprintf('===============================================================\n');
fprintf('       PIXEL-WISE NRMSE RUN-TO-RUN VARIABILITY STATISTICS\n');
fprintf('===============================================================\n');
fprintf('%-8s %-14s %-14s %-14s %-14s\n', ...
    'Model', ...
    'Mean STD', ...
    'Median STD', ...
    'P10 STD', ...
    'P90 STD');
fprintf('---------------------------------------------------------------\n');

for iModel = 1:nModel
    fprintf('%-8s %-14.6f %-14.6f %-14.6f %-14.6f\n', ...
        modelNames{iModel}, ...
        meanSTD(iModel), ...
        medianSTD(iModel), ...
        percentile10STD(iModel), ...
        percentile90STD(iModel));
end

fprintf('===============================================================\n');

% Also display the MATLAB table
disp(uncertaintyTable);

% ------------------- Plot Fig1-Fig6 -------------------
for iModel = 1:nModel
    modelName = modelNames{iModel};
    stdMap = NRMSE_STD_maps(:, :, iModel);

    figureTitle = sprintf('Global NRMSE-STD of %s model', modelName);

    figureFile = fullfile(outDir, ...
        sprintf('Fig%d_Global_NRMSE_STD_%s_model.png', iModel, modelName));

    plotGlobalSTDMap( ...
        iModel, ...
        Lon2D, ...
        Lat2D, ...
        stdMap, ...
        figureTitle, ...
        figureFile, ...
        figPosition, ...
        uncertaintyColormap, ...
        fontName, ...
        titleFontSize, ...
        axisLabelFontSize, ...
        tickFontSize, ...
        colorbarFontSize, ...
        colorbarLabelFontSize, ...
        outputResolution);
end

% ------------------- Save statistics table -------------------
csvFile = fullfile(outDir, ...
    's1_overall_performance3_NRMSE_STD_statistics.csv');

writetable(uncertaintyTable, csvFile);

% ------------------- Save calculation results -------------------
matFile = fullfile(outDir, 's1_overall_performance3_results.mat');

save(matFile, ...
    'modelNames', ...
    'modelFiles', ...
    'nrmseVarNames', ...
    'NRMSE_STD_maps', ...
    'meanSTD', ...
    'medianSTD', ...
    'percentile10STD', ...
    'percentile90STD', ...
    'validPixelCounts', ...
    'uncertaintyTable', ...
    'validMask', ...
    'lat', ...
    'lon', ...
    'figPosition');

fprintf('\nCompleted: pixel-wise NRMSE uncertainty analysis.\n');
fprintf('Output directory:\n%s\n', outDir);

%
% ==================== 4. Additional Functions ====================

function plotGlobalSTDMap( ...
    figureNumber, ...
    Lon2D, ...
    Lat2D, ...
    stdMap, ...
    figureTitle, ...
    figureFile, ...
    figPosition, ...
    uncertaintyColormap, ...
    fontName, ...
    titleFontSize, ...
    axisLabelFontSize, ...
    tickFontSize, ...
    colorbarFontSize, ...
    colorbarLabelFontSize, ...
    outputResolution)

    %% ---------- Create figure ----------

    fig = figure(figureNumber);
    clf(fig);

    set(fig, ...
        'Color', 'w', ...
        'Position', figPosition);

    % Explicitly create an axes object
    ax = axes('Parent', fig);

    %% ---------- Plot global NRMSE standard-deviation map ----------

    pcolor(ax, Lon2D, Lat2D, stdMap);
    shading(ax, 'flat');

    %% ---------- Axis settings ----------

    set(ax, ...
        'YDir', 'normal', ...
        'FontName', fontName, ...
        'FontSize', tickFontSize, ...
        'LineWidth', 0.8);

    axis(ax, 'equal');
    axis(ax, 'tight');

    xlim(ax, [min(Lon2D(:)), max(Lon2D(:))]);
    ylim(ax, [min(Lat2D(:)), max(Lat2D(:))]);

    xlabel(ax, 'Longitude', ...
        'FontName', fontName, ...
        'FontSize', axisLabelFontSize);

    ylabel(ax, 'Latitude', ...
        'FontName', fontName, ...
        'FontSize', axisLabelFontSize);

    title(ax, figureTitle, ...
        'FontName', fontName, ...
        'FontSize', titleFontSize, ...
        'FontWeight', 'normal');

    %% ---------- Colormap and color limits ----------

    colormap(ax, uncertaintyColormap);

    % A shared fixed color range makes all six maps directly comparable.
    clim(ax, [0, 0.04]);

    %% ---------- Colorbar settings ----------

    cb = colorbar(ax);
    cb.FontName = fontName;
    cb.FontSize = colorbarFontSize;

    ylabel(cb, 'NRMSE standard deviation', ...
        'FontName', fontName, ...
        'FontSize', colorbarLabelFontSize);

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', outputResolution);
end
