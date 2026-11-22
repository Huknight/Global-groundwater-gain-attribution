%% SARIMA - Global univariate groundwater forecasting
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Forecast monthly groundwater-storage change (dGWS) at each valid
%   1-degree grid cell using a seasonal autoregressive integrated moving-average (SARIMA) model.
%
% Input
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   dGWS_valid: [grid, month], January 2004 to December 2024.
%
% Rolling evaluation
%   Test 2022: train 2004-2020, validate 2021.
%   Test 2023: train 2004-2021, validate 2022.
%   Test 2024: train 2004-2022, validate 2023.
%   Each forecast covers 12 months from the preceding December.
%
% Model
%   SARIMA(p,0,1)(1,1,1)_12, p = [1, 6, 12, 24, 36, 48, 60].
%   The 12-month seasonal difference is set by Seasonality = 12.
%   Seasonal AR and MA terms are specified at lag 12.
%   MATLAB estimate fits unknown parameters by Gaussian maximum likelihood.
%   The validation year selects p using NRMSE. Model coefficients are
%   estimated from the training period only and are not refitted after
%   validation. Observed training and validation history initializes the
%   AR and MA terms at the test forecast origin.
%
% Evaluation
%   NRMSE = RMSE / (max(test dGWS) - min(test dGWS)).
%   TDA compares directions of consecutive changes in monthly dGWS values:
%   sign(dGWS(t) - dGWS(t-1)), with 11 comparisons per test year.
%   TDA is reported in percent and is not used for model selection.
%
% Output
%   SARIMA_Pred:                 [lat, lon, 36], test years 2022-2024
%   SARIMA_Val_Pred:             [lat, lon, 36], validation years 2021-2023
%   SARIMA_NRMSE, SARIMA_TDA:        [lat, lon, 3], one layer per test year
%   SARIMA_best_p:               [lat, lon, 3], validation-selected nonseasonal AR order
%   SARIMA_Val_NRMSE_best:       [lat, lon, 3]
%   SARIMA_TestSuccess:          [lat, lon, 3], successful forecasts
%   SARIMA_FailureCode:          [lat, lon, 3], see Meta.failureCodes
%   Outputs also include valid-grid matrices, coordinates and metadata.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
clearvars;
clc;

projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
dataFile = fullfile(projectRoot, 's1_master_database', '03_master_database', ...
    'GLDAS_common_valid_mask_1deg_200401_202412.mat');
outFolder = fullfile(projectRoot, 's1_master_database', '04_linear_models');
outFileName = 'SARIMA_univariate_global_200401_202412.mat';
checkpointFileName = 'SARIMA_univariate_global_checkpoint.mat';

% Existing final output: "overwrite" or "stop".
outputMode = "overwrite";

% Resume only when the checkpoint matches the input data and model settings.
resumeFromCheckpoint = true;

useParallel = true;
numWorkers = 15
















;
block_num = 10;

p_list = [1, 6, 12, 24, 36, 48, 60];
testYears = [2022, 2023, 2024];
forecastHorizon = 12;
d = 0;
q = 1;
S = 12;
P = 1;
D_seasonal = 1;
Q = 1;

makeDiagnosticPlot = true;
plotNRMSELimits = [0, 0.4];
checkpointVersion = 1;

%% =========================================================================
% 3. MAIN EXECUTION
% =========================================================================

%% 3.1 Check settings and read the master database
if ~isfile(dataFile)
    error('Master database not found: %s', dataFile);
end
if ~ismember(outputMode, ["overwrite", "stop"])
    error('outputMode must be "overwrite" or "stop".');
end
if ~isscalar(block_num) || block_num < 1 || block_num ~= floor(block_num)
    error('block_num must be a positive integer.');
end
if useParallel && (~isscalar(numWorkers) || numWorkers < 1 || ...
        numWorkers ~= floor(numWorkers))
    error('numWorkers must be a positive integer.');
end
if any(p_list < 1) || any(p_list ~= floor(p_list)) || ...
        numel(unique(p_list)) ~= numel(p_list)
    error('p_list must contain distinct positive integers.');
end
if ~isscalar(d) || ~isscalar(q) || d ~= 0 || q ~= 1 || ...
        S ~= 12 || P ~= 1 || D_seasonal ~= 1 || Q ~= 1
    error('This experiment requires SARIMA(p,0,1)(1,1,1)_12.');
end
if ~isequal(testYears, [2022, 2023, 2024]) || forecastHorizon ~= 12
    error('This experiment requires three annual 12-month forecasts, 2022-2024.');
end

if ~isfolder(outFolder)
    [ok, msg] = mkdir(outFolder);
    if ~ok
        error('Cannot create output directory: %s\n%s', outFolder, msg);
    end
end
outFile = fullfile(outFolder, outFileName);
checkpointFile = fullfile(outFolder, checkpointFileName);
if isfile(outFile) && outputMode == "stop"
    error('Final output already exists: %s', outFile);
end

requiredVars = {'dGWS_valid', 'CommonMask', 'CommonValid', ...
    'validLinearIndex', 'validLat', 'validLon', 'lat', 'lon', 'ym'};
checkVarsExist(who('-file', dataFile), requiredVars);

D = load(dataFile, requiredVars{:});
dGWS_valid = double(D.dGWS_valid);
CommonMask = D.CommonMask;
CommonValid = logical(D.CommonValid);
validLinearIndex = double(D.validLinearIndex(:));
validLat = double(D.validLat(:));
validLon = double(D.validLon(:));
lat = double(D.lat(:));
lon = double(D.lon(:));
ym = double(D.ym(:));
clear D

nLat = numel(lat);
nLon = numel(lon);
nValid = size(dGWS_valid, 1);
nTime = numel(ym);
expectedYM = makeYMList(2004, 1, 2024, 12);

if nLat ~= 150 || nLon ~= 360 || ~isequal(ym, expectedYM)
    error('Expected a 150 x 360 grid and monthly data from 2004-01 to 2024-12.');
end
if ~isequal(size(dGWS_valid), [nValid, nTime]) || nTime ~= 252 || ...
        any(~isfinite(dGWS_valid(:)))
    error('dGWS_valid must contain finite [grid, 252] time series.');
end
if ~isequal(size(CommonValid), [nLat, nLon]) || ...
        ~isequal(size(CommonMask), [nLat, nLon]) || ...
        nnz(CommonValid) ~= nValid || ...
        ~isequal(validLinearIndex, find(CommonValid))
    error('The common mask does not match the valid-grid time series.');
end
if any(CommonMask(CommonValid) ~= 1) || ...
        any(~isnan(CommonMask(~CommonValid))) || ...
        numel(validLat) ~= nValid || numel(validLon) ~= nValid
    error('CommonMask or valid-grid coordinates are inconsistent.');
end
[Lon2D, Lat2D] = meshgrid(lon, lat);
if any(abs(Lat2D(validLinearIndex) - validLat) > 1e-8) || ...
        any(abs(Lon2D(validLinearIndex) - validLon) > 1e-8)
    error('Valid-grid coordinate vectors do not match validLinearIndex.');
end

nRound = numel(testYears);
nPredMonth = nRound * forecastHorizon;
predYM = makeYMList(testYears(1), 1, testYears(end), 12);
valYM = makeYMList(testYears(1)-1, 1, testYears(end)-1, 12);
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);

block_size = max(1, ceil(nValid / block_num));
total_blocks = ceil(nValid / block_size);
fprintf('SARIMA model: %d grid cells, %d months, %d test rounds.\n', ...
    nValid, nTime, nRound);
fprintf('SARIMA(p,%d,%d)(%d,%d,%d)_%d, p = %s\n', ...
    d, q, P, D_seasonal, Q, S, mat2str(p_list));
fprintf('Blocks: %d; approximately %d grid cells per block.\n', ...
    total_blocks, block_size);

%% 3.2 Prepare checkpoint and result arrays
% SHA-256 protects against resuming results from another input or setup.
runSignature = computeRunSignature(dGWS_valid, validLinearIndex, ...
    lat, lon, ym, p_list, d, q, S, P, D_seasonal, Q, ...
    testYears, forecastHorizon, block_num, checkpointVersion);

if isfile(checkpointFile) && ~resumeFromCheckpoint
    delete(checkpointFile);
end

if resumeFromCheckpoint && isfile(checkpointFile)
    fprintf('Reading checkpoint: %s\n', checkpointFile);
    C = load(checkpointFile);
    validateCheckpoint(C, runSignature, nValid, nRound, nPredMonth, total_blocks);

    SARIMA_Pred_valid = C.SARIMA_Pred_valid;
    SARIMA_Val_Pred_valid = C.SARIMA_Val_Pred_valid;
    SARIMA_NRMSE_valid = C.SARIMA_NRMSE_valid;
    SARIMA_TDA_valid = C.SARIMA_TDA_valid;
    SARIMA_best_p_valid = C.SARIMA_best_p_valid;
    SARIMA_Val_NRMSE_best_valid = C.SARIMA_Val_NRMSE_best_valid;
    SARIMA_RMSE_valid = C.SARIMA_RMSE_valid;
    SARIMA_TestRange_valid = C.SARIMA_TestRange_valid;
    SARIMA_TestSuccess_valid = C.SARIMA_TestSuccess_valid;
    SARIMA_FailureCode_valid = C.SARIMA_FailureCode_valid;
    SARIMA_FailureID_valid = C.SARIMA_FailureID_valid;
    blockDone = C.blockDone;
    cumulativeComputeTimeSeconds = C.cumulativeComputeTimeSeconds;
    clear C
else
    SARIMA_Pred_valid = NaN(nValid, nPredMonth);
    SARIMA_Val_Pred_valid = NaN(nValid, nPredMonth);
    SARIMA_NRMSE_valid = NaN(nValid, nRound);
    SARIMA_TDA_valid = NaN(nValid, nRound);
    SARIMA_best_p_valid = NaN(nValid, nRound);
    SARIMA_Val_NRMSE_best_valid = NaN(nValid, nRound);
    SARIMA_RMSE_valid = NaN(nValid, nRound);
    SARIMA_TestRange_valid = NaN(nValid, nRound);
    SARIMA_TestSuccess_valid = false(nValid, nRound);
    SARIMA_FailureCode_valid = uint8(zeros(nValid, nRound));
    SARIMA_FailureID_valid = cell(nValid, nRound);
    blockDone = false(total_blocks, 1);
    cumulativeComputeTimeSeconds = 0;
end

%% 3.3 Start parallel pool
if useParallel
    try
        poolObj = gcp('nocreate');
        if isempty(poolObj)
            parpool('local', numWorkers);
        elseif poolObj.NumWorkers ~= numWorkers
            delete(poolObj);
            parpool('local', numWorkers);
        end
    catch ME
        warning('SKLGB:ParallelPoolUnavailable', ...
            'Parallel pool unavailable; using serial processing: %s', ME.message);
        useParallel = false;
    end
end

if useParallel
    Y_valid_const = parallel.pool.Constant(dGWS_valid);
end

%% 3.4 Run global SARIMA forecasts
fprintf('\nStarting SARIMA computation...\n');
totalTic = tic;

for block = 1:total_blocks
    if blockDone(block)
        fprintf('Block %d/%d already complete.\n', block, total_blocks);
        continue;
    end

    start_idx = (block - 1) * block_size + 1;
    end_idx = min(block * block_size, nValid);
    currentRows = start_idx:end_idx;
    block_len = numel(currentRows);

    blockPred = NaN(block_len, nPredMonth);
    blockValPred = NaN(block_len, nPredMonth);
    blockNRMSE = NaN(block_len, nRound);
    blockTDA = NaN(block_len, nRound);
    blockBestP = NaN(block_len, nRound);
    blockValBest = NaN(block_len, nRound);
    blockRMSE = NaN(block_len, nRound);
    blockTestRange = NaN(block_len, nRound);
    blockSuccess = false(block_len, nRound);
    blockFailureCode = uint8(zeros(block_len, nRound));
    blockFailureID = cell(block_len, nRound);

    blockTic = tic;
    if useParallel
        parfor i = 1:block_len
            Y_pixel = Y_valid_const.Value(currentRows(i), :)';
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ok3, code3, id3] = ...
                runSARIMAOnePixelAllRounds(Y_pixel, p_list, d, q, ...
                S, P, D_seasonal, Q, testYears, forecastHorizon);
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockValBest(i, :) = vn3(:)';
            blockRMSE(i, :) = rmse3(:)';
            blockTestRange(i, :) = range3(:)';
            blockSuccess(i, :) = ok3(:)';
            blockFailureCode(i, :) = code3(:)';
            blockFailureID(i, :) = id3(:)';
        end
    else
        for i = 1:block_len
            Y_pixel = dGWS_valid(currentRows(i), :)';
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ok3, code3, id3] = ...
                runSARIMAOnePixelAllRounds(Y_pixel, p_list, d, q, ...
                S, P, D_seasonal, Q, testYears, forecastHorizon);
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockValBest(i, :) = vn3(:)';
            blockRMSE(i, :) = rmse3(:)';
            blockTestRange(i, :) = range3(:)';
            blockSuccess(i, :) = ok3(:)';
            blockFailureCode(i, :) = code3(:)';
            blockFailureID(i, :) = id3(:)';
        end
    end

    blockElapsed = toc(blockTic);
    cumulativeComputeTimeSeconds = cumulativeComputeTimeSeconds + blockElapsed;

    SARIMA_Pred_valid(currentRows, :) = blockPred;
    SARIMA_Val_Pred_valid(currentRows, :) = blockValPred;
    SARIMA_NRMSE_valid(currentRows, :) = blockNRMSE;
    SARIMA_TDA_valid(currentRows, :) = blockTDA;
    SARIMA_best_p_valid(currentRows, :) = blockBestP;
    SARIMA_Val_NRMSE_best_valid(currentRows, :) = blockValBest;
    SARIMA_RMSE_valid(currentRows, :) = blockRMSE;
    SARIMA_TestRange_valid(currentRows, :) = blockTestRange;
    SARIMA_TestSuccess_valid(currentRows, :) = blockSuccess;
    SARIMA_FailureCode_valid(currentRows, :) = blockFailureCode;
    SARIMA_FailureID_valid(currentRows, :) = blockFailureID;
    blockDone(block) = true;

    checkpointInfo = struct('version', checkpointVersion, ...
        'runSignature', runSignature);
    checkpointTemp = fullfile(outFolder, [checkpointFileName, '.tmp.mat']);
    save(checkpointTemp, ...
        'SARIMA_Pred_valid', 'SARIMA_Val_Pred_valid', ...
        'SARIMA_NRMSE_valid', 'SARIMA_TDA_valid', ...
        'SARIMA_best_p_valid', 'SARIMA_Val_NRMSE_best_valid', ...
        'SARIMA_RMSE_valid', 'SARIMA_TestRange_valid', ...
        'SARIMA_TestSuccess_valid', 'SARIMA_FailureCode_valid', 'SARIMA_FailureID_valid', ...
        'blockDone', 'cumulativeComputeTimeSeconds', 'checkpointInfo', '-v7.3');
    [moved, msg] = movefile(checkpointTemp, checkpointFile, 'f');
    if ~moved
        error('Checkpoint replacement failed: %s', msg);
    end

    completedBlocks = nnz(blockDone);
    elapsedThisRun = toc(totalTic);
    fprintf('Block %d/%d | %.2f min | %d/%d completed | %.1f%% | run %.2f min\n', ...
        block, total_blocks, blockElapsed / 60, completedBlocks, total_blocks, ...
        100 * completedBlocks / total_blocks, elapsedThisRun / 60);
end

if ~all(blockDone)
    error('Not all blocks were completed.');
end
elapsedComputeThisRun = toc(totalTic);
fprintf('\nSARIMA computation finished.\n');
fprintf('Computation time this run: %.2f minutes (%.2f hours)\n', ...
    elapsedComputeThisRun / 60, elapsedComputeThisRun / 3600);
fprintf('Cumulative block computation time: %.2f minutes (%.2f hours)\n', ...
    cumulativeComputeTimeSeconds / 60, cumulativeComputeTimeSeconds / 3600);
if useParallel
    delete(Y_valid_const);
end

%% 3.5 Map valid-grid results to 1-degree fields
SARIMA_Pred = mapValidToGrid3D(SARIMA_Pred_valid, validLinearIndex, nLat, nLon);
SARIMA_Val_Pred = mapValidToGrid3D(SARIMA_Val_Pred_valid, validLinearIndex, nLat, nLon);
SARIMA_NRMSE = mapValidToGrid3D(SARIMA_NRMSE_valid, validLinearIndex, nLat, nLon);
SARIMA_TDA = mapValidToGrid3D(SARIMA_TDA_valid, validLinearIndex, nLat, nLon);
SARIMA_best_p = mapValidToGrid3D(SARIMA_best_p_valid, validLinearIndex, nLat, nLon);
SARIMA_Val_NRMSE_best = mapValidToGrid3D(SARIMA_Val_NRMSE_best_valid, validLinearIndex, nLat, nLon);
SARIMA_RMSE = mapValidToGrid3D(SARIMA_RMSE_valid, validLinearIndex, nLat, nLon);
SARIMA_TestRange = mapValidToGrid3D(SARIMA_TestRange_valid, validLinearIndex, nLat, nLon);
SARIMA_TestSuccess = mapValidToGrid3D(double(SARIMA_TestSuccess_valid), ...
    validLinearIndex, nLat, nLon);
SARIMA_FailureCode = mapValidToGrid3D(double(SARIMA_FailureCode_valid), ...
    validLinearIndex, nLat, nLon);

SARIMA_NRMSE_mean = mean(SARIMA_NRMSE, 3, 'omitnan');
SARIMA_TDA_mean = mean(SARIMA_TDA, 3, 'omitnan');
SARIMA_Val_NRMSE_best_mean = mean(SARIMA_Val_NRMSE_best, 3, 'omitnan');

% CompleteThreeYears identifies cells available for strict model comparison.
CompleteThreeYears_valid = all(SARIMA_TestSuccess_valid, 2);
SARIMA_CompleteThreeYears = false(nLat, nLon);
SARIMA_CompleteThreeYears(validLinearIndex) = CompleteThreeYears_valid;
SARIMA_NRMSE_mean_complete = SARIMA_NRMSE_mean;
SARIMA_NRMSE_mean_complete(~SARIMA_CompleteThreeYears) = NaN;

fprintf('\nSARIMA global results:\n');
fprintf('Completed test forecasts: %d / %d grid-year cases.\n', ...
    nnz(SARIMA_TestSuccess_valid), numel(SARIMA_TestSuccess_valid));
fprintf('Complete three-year grid cells: %d / %d.\n', ...
    nnz(CompleteThreeYears_valid), nValid);
for r = 1:nRound
    fprintf('  %d success: %d / %d\n', testYears(r), ...
        nnz(SARIMA_TestSuccess_valid(:, r)), nValid);
    for code = 1:4
        nFailure = nnz(SARIMA_FailureCode_valid(:, r) == code);
        if nFailure > 0
            fprintf('    Failure code %d: %d\n', code, nFailure);
        end
    end
end
fprintf('SARIMA_NRMSE_mean (finite values): %.4f\n', ...
    mean(SARIMA_NRMSE_mean(CommonValid), 'omitnan'));
fprintf('SARIMA_TDA_mean (finite values): %.2f%%\n', ...
    mean(SARIMA_TDA_mean(CommonValid), 'omitnan'));

%% 3.6 Save output
Meta = struct();
Meta.model = 'SARIMA';
Meta.input = 'Monthly dGWS from the common GLDAS valid-grid database';
Meta.period = '2004-01 to 2024-12';
Meta.predictionPeriod = '2022-01 to 2024-12';
Meta.validationPeriod = '2021-01 to 2023-12';
Meta.arrayConvention = 'Maps: [latitude, longitude, time]; valid series: [grid, time]';
Meta.rollingStrategy = ['Train 2004 through test year minus two; ' ...
    'validate in the following year; test in the subsequent year'];
Meta.sarimaEstimation = 'MATLAB arima/estimate, Gaussian maximum likelihood';
Meta.p_list = p_list;
Meta.d = d;
Meta.q = q;
Meta.S = S;
Meta.P = P;
Meta.D_seasonal = D_seasonal;
Meta.Q = Q;
Meta.modelOrder = 'SARIMA(p,0,1)(1,1,1)_12';
Meta.finalPrediction = ['No refit using validation data; selected SARIMA parameters ' ...
    'are estimated on training data; observed validation history ' ...
    'is available at the test forecast origin'];
Meta.NRMSE = 'Annual RMSE divided by observed test-year dGWS range';
Meta.TDA = ['Percentage of matching signs of successive month-to-month ' ...
    'differences in dGWS within a test year (11 comparisons)'];
Meta.failureCodes = ['0 = success; 1 = no finite validation candidate; ' ...
    '2 = test forecasting error; 3 = invalid test predictions; ' ...
    '4 = nonfinite test NRMSE'];
Meta.SARIMA_FailureID = 'MATLAB exception identifiers for failed forecast attempts';
Meta.meanMetricRule = 'Three-year mean computed with omitnan';
Meta.completeMetricRule = ['SARIMA_NRMSE_mean_complete requires successful forecasts ' ...
    'in all three test years'];
Meta.nValid = nValid;
Meta.block_num = block_num;
Meta.numWorkers = numWorkers;
Meta.usedParallel = useParallel;
Meta.runSignature = runSignature;
Meta.computeTimeMinutes = cumulativeComputeTimeSeconds / 60;
Meta.computeTimeThisRunMinutes = elapsedComputeThisRun / 60;
Meta.scriptName = 'SARIMA_global_forecast.m';

finalTemp = fullfile(outFolder, [outFileName, '.tmp.mat']);
fprintf('\nSaving final MAT-file: %s\n', outFile);
save(finalTemp, ...
    'SARIMA_Pred', 'SARIMA_Val_Pred', ...
    'SARIMA_NRMSE', 'SARIMA_NRMSE_mean', 'SARIMA_NRMSE_mean_complete', ...
    'SARIMA_TDA', 'SARIMA_TDA_mean', ...
    'SARIMA_best_p', 'SARIMA_Val_NRMSE_best', 'SARIMA_Val_NRMSE_best_mean', ...
    'SARIMA_RMSE', 'SARIMA_TestRange', ...
    'SARIMA_TestSuccess', 'SARIMA_FailureCode', 'SARIMA_CompleteThreeYears', ...
    'SARIMA_Pred_valid', 'SARIMA_Val_Pred_valid', ...
    'SARIMA_NRMSE_valid', 'SARIMA_TDA_valid', 'SARIMA_best_p_valid', ...
    'SARIMA_Val_NRMSE_best_valid', 'SARIMA_RMSE_valid', 'SARIMA_TestRange_valid', ...
    'SARIMA_TestSuccess_valid', 'SARIMA_FailureCode_valid', 'SARIMA_FailureID_valid', ...
    'p_list', 'd', 'q', 'S', 'P', 'D_seasonal', 'Q', ...
    'testYears', 'predYM', 'valYM', ...
    'CommonMask', 'CommonValid', 'validLinearIndex', 'nValid', ...
    'validLat', 'validLon', 'lat', 'lon', 'ym', 'yearVec', 'monthVec', ...
    'Meta', '-v7.3');
[moved, msg] = movefile(finalTemp, outFile, 'f');
if ~moved
    error('Final MAT-file replacement failed: %s', msg);
end
fprintf('Saved successfully.\n');

%% 3.7 Diagnostic map
if makeDiagnosticPlot
    figure('Name', 'SARIMA mean NRMSE');
    p = imagesc(lon, lat, SARIMA_NRMSE_mean);
    set(p, 'AlphaData', double(isfinite(SARIMA_NRMSE_mean)));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(jet);
    clim(plotNRMSELimits);
    cb = colorbar;
    ylabel(cb, 'NRMSE');
    title('SARIMA mean NRMSE, test years 2022-2024');
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

elapsedTotalSeconds = toc(totalTic);
fprintf('\nSARIMA processing completed.\n');
fprintf('Total time since computation started (including output): %.2f minutes (%.2f hours)\n', ...
    elapsedTotalSeconds / 60, elapsedTotalSeconds / 3600);

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

function [pred36, valPred36, nrmse3, tda3, bestP3, valBest3, ...
    rmse3, testRange3, success3, failureCode3, failureID3] = ...
    runSARIMAOnePixelAllRounds(Y, p_list, d, q, S, P, D_seasonal, ...
    Q, testYears, horizon)

    Y = double(Y(:));
    nRound = numel(testYears);
    pred36 = NaN(nRound * horizon, 1);
    valPred36 = NaN(nRound * horizon, 1);
    nrmse3 = NaN(nRound, 1);
    tda3 = NaN(nRound, 1);
    bestP3 = NaN(nRound, 1);
    valBest3 = NaN(nRound, 1);
    rmse3 = NaN(nRound, 1);
    testRange3 = NaN(nRound, 1);
    success3 = false(nRound, 1);
    failureCode3 = uint8(zeros(nRound, 1));
    failureID3 = cell(nRound, 1);

    for r = 1:nRound
        testYear = testYears(r);
        trainEndIdx = (testYear - 2005) * 12;
        valIdx = (trainEndIdx + 1):(trainEndIdx + horizon);
        testIdx = (trainEndIdx + horizon + 1): ...
            (trainEndIdx + 2 * horizon);
        predIdx = ((r - 1) * horizon + 1):(r * horizon);

        Y_train = Y(1:trainEndIdx);
        Y_val = Y(valIdx);
        Y_test = Y(testIdx);
        bestValNRMSE = Inf;
        bestModel = [];
        bestValPred = [];
        bestP = NaN;
        firstValidationError = '';

        for k = 1:numel(p_list)
            p = p_list(k);
            maxLag = max([p, q, S * P, S * Q]);
            if maxLag >= numel(Y_train)
                continue;
            end
            try
                model = makeSARIMAModel(p, d, q, S, P, D_seasonal, Q);
                estModel = estimate(model, Y_train, 'Display', 'off');
                Y_val_pred = forecastSARIMANumeric(estModel, Y_train, horizon);
                valScore = calcNRMSE(Y_val, Y_val_pred);
                if isfinite(valScore) && valScore < bestValNRMSE
                    bestValNRMSE = valScore;
                    bestModel = estModel;
                    bestValPred = Y_val_pred;
                    bestP = p;
                end
            catch ME
                if isempty(firstValidationError)
                    firstValidationError = ME.identifier;
                end
            end
        end

        if isempty(bestModel)
            failureCode3(r) = 1;
            failureID3{r} = firstValidationError;
            continue;
        end

        bestP3(r) = bestP;
        valBest3(r) = bestValNRMSE;
        valPred36(predIdx) = bestValPred;

        % Use the observed validation year only as forecast history.
        Y_train_val = Y(1:trainEndIdx + horizon);
        try
            Y_test_pred = forecastSARIMANumeric(bestModel, Y_train_val, horizon);
        catch ME
            failureCode3(r) = 2;
            failureID3{r} = ME.identifier;
            continue;
        end

        if any(~isfinite(Y_test_pred))
            failureCode3(r) = 3;
            continue;
        end

        pred36(predIdx) = Y_test_pred;
        testNRMSE = calcNRMSE(Y_test, Y_test_pred);
        if ~isfinite(testNRMSE)
            failureCode3(r) = 4;
            continue;
        end

        nrmse3(r) = testNRMSE;
        rmse3(r) = sqrt(mean((Y_test_pred - Y_test).^2));
        testRange3(r) = max(Y_test) - min(Y_test);
        tda3(r) = calcTDA(Y_test, Y_test_pred);
        success3(r) = true;
    end
end

function model = makeSARIMAModel(p, d, q, S, P, D_seasonal, Q)
    if D_seasonal == 1
        model = arima('ARLags', 1:p, 'D', d, 'MALags', 1:q, ...
            'Seasonality', S, 'SARLags', S * (1:P), ...
            'SMALags', S * (1:Q));
    else
        model = arima('ARLags', 1:p, 'D', d, 'MALags', 1:q, ...
            'SARLags', S * (1:P), 'SMALags', S * (1:Q));
    end
end

function yhat = forecastSARIMANumeric(model, historyY, horizon)
    yhat = forecast(model, horizon, 'Y0', historyY(:));
    yhat = double(yhat(:));
    if numel(yhat) ~= horizon
        error('SKLGB:ForecastLength', ...
            'Forecast returned %d months; expected %d.', numel(yhat), horizon);
    end
end

function nrmse = calcNRMSE(obs, pred)
    obs = double(obs(:));
    pred = double(pred(:));
    if isempty(obs) || numel(obs) ~= numel(pred) || ...
            any(~isfinite(obs)) || any(~isfinite(pred))
        nrmse = NaN;
        return;
    end
    rmse = sqrt(mean((pred - obs).^2));
    denom = max(obs) - min(obs);
    if denom == 0
        if rmse == 0
            nrmse = 0;
        else
            nrmse = Inf;
        end
    else
        nrmse = rmse / denom;
    end
end

function tda = calcTDA(obs, pred)
    obs = double(obs(:));
    pred = double(pred(:));
    if numel(obs) ~= numel(pred) || numel(obs) < 2 || ...
            any(~isfinite(obs)) || any(~isfinite(pred))
        tda = NaN;
        return;
    end
    actualDiff = diff(obs);
    predDiff = diff(pred);
    tda = mean(sign(actualDiff) == sign(predDiff)) * 100;
end

function Ygrid = mapValidToGrid3D(Yvalid, validLinearIndex, nLat, nLon)
    [~, nLayer] = size(Yvalid);
    Ygrid = NaN(nLat, nLon, nLayer);
    for k = 1:nLayer
        tmp = NaN(nLat, nLon);
        tmp(validLinearIndex) = Yvalid(:, k);
        Ygrid(:, :, k) = tmp;
    end
end

function ym = makeYMList(y1, m1, y2, m2)
    d1 = datetime(y1, m1, 1);
    d2 = datetime(y2, m2, 1);
    if d2 < d1
        error('End date precedes start date.');
    end
    d = d1:calmonths(1):d2;
    ym = double(year(d(:)) * 100 + month(d(:)));
end

function checkVarsExist(existing, required)
    for i = 1:numel(required)
        if ~any(strcmp(existing, required{i}))
            error('Missing variable in master database: %s', required{i});
        end
    end
end

function signature = computeRunSignature(data, index, lat, lon, ym, ...
    p_list, d, q, S, P, D_seasonal, Q, ...
    testYears, horizon, block_num, version)
    md = java.security.MessageDigest.getInstance('SHA-256');
    updateDigest(md, 'SARIMA_global_forecast');
    updateDigest(md, double(version));
    updateDigest(md, data);
    updateDigest(md, index);
    updateDigest(md, lat);
    updateDigest(md, lon);
    updateDigest(md, ym);
    updateDigest(md, double(p_list));
    updateDigest(md, double(d));
    updateDigest(md, double(q));
    updateDigest(md, double(S));
    updateDigest(md, double(P));
    updateDigest(md, double(D_seasonal));
    updateDigest(md, double(Q));
    updateDigest(md, double(testYears));
    updateDigest(md, double(horizon));
    updateDigest(md, double(block_num));
    hashBytes = typecast(int8(md.digest()), 'uint8');
    signature = lower(reshape(dec2hex(hashBytes, 2).', 1, []));
end

function updateDigest(md, value)
    % Include the array type and shape to make the signature unambiguous.
    if ischar(value) || isstring(value)
        bytes = int8(uint8(char(value)));
        dims = uint64(numel(bytes));
        md.update(int8(uint8('text')));
        md.update(typecast(dims, 'int8'));
        md.update(bytes);
    else
        v = double(value);
        md.update(int8(uint8('double')));
        md.update(typecast(uint64(size(v)), 'int8'));
        md.update(typecast(v(:), 'int8'));
    end
end

function validateCheckpoint(C, signature, nValid, nRound, nPredMonth, totalBlocks)
    required = {'checkpointInfo', 'blockDone', ...
        'SARIMA_Pred_valid', 'SARIMA_Val_Pred_valid', ...
        'SARIMA_NRMSE_valid', 'SARIMA_TDA_valid', ...
        'SARIMA_best_p_valid', 'SARIMA_Val_NRMSE_best_valid', ...
        'SARIMA_RMSE_valid', 'SARIMA_TestRange_valid', ...
        'SARIMA_TestSuccess_valid', 'SARIMA_FailureCode_valid', ...
        'SARIMA_FailureID_valid', 'cumulativeComputeTimeSeconds'};
    for k = 1:numel(required)
        if ~isfield(C, required{k})
            error('Checkpoint is missing %s. Use resumeFromCheckpoint=false.', ...
                required{k});
        end
    end
    if ~isfield(C.checkpointInfo, 'runSignature') || ...
            ~strcmp(C.checkpointInfo.runSignature, signature)
        error(['Checkpoint input data or model settings differ from this run. ' ...
            'Use resumeFromCheckpoint=false to start a new computation.']);
    end
    matrix3 = {'SARIMA_NRMSE_valid', 'SARIMA_TDA_valid', 'SARIMA_best_p_valid', ...
        'SARIMA_Val_NRMSE_best_valid', 'SARIMA_RMSE_valid', 'SARIMA_TestRange_valid', ...
        'SARIMA_TestSuccess_valid', 'SARIMA_FailureCode_valid', 'SARIMA_FailureID_valid'};
    for k = 1:numel(matrix3)
        if ~isequal(size(C.(matrix3{k})), [nValid, nRound])
            error('Checkpoint contains an invalid %s array.', matrix3{k});
        end
    end
    if ~isequal(size(C.SARIMA_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.SARIMA_Val_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.blockDone), [totalBlocks, 1])
        error('Checkpoint dimensions do not match the current run.');
    end
end
