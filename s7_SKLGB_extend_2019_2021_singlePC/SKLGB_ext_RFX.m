function SKLGB_ext_RFX(cfg)
% PILOT WRAPPER: adapted from RFX_global_forecast.m.
% Only file paths, pilot seed IDs, run settings and the two tested
% hyperparameters are changed; the original model code follows.
% Each configuration searches the FULL original validation p/x grid.

%% RFX - Global multivariate groundwater forecasting
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Forecast monthly groundwater-storage change (dGWS) at each valid
%   1-degree grid cell using bagged regression trees and historical
%   hydroclimatic predictors.
%
% Input
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   dGWS, precipitation, evaporation, runoff, soil-moisture change and
%   surface temperature: [grid, month], January 2004 to December 2024.
%
% Rolling evaluation
%   Test 2022: train 2004-2020, validate 2021.
%   Test 2023: train 2004-2021, validate 2022.
%   Test 2024: train 2004-2022, validate 2023.
%   The forecast origin precedes each 12-month test year.
%
% Model
%   Bagged regression trees, fitrensemble with Method = 'Bag'.
%   50 trees, MinLeafSize = 5; no predictor or target standardization.
%   dGWS history p: [1, 6, 12, 24, 36, 48, 60] months.
%   Exogenous history x: [1, 3, 6] months for each of five predictors.
%   NumVariablesToSample = max(1, floor((p + 5*x) / 3)).
%   For each (p, x) and replicate, twelve separate ensembles are trained
%   to predict horizons h = 1,...,12 from a common forecast origin.
%   Validation-year NRMSE selects p and x; selected fitted ensembles are
%   reused for test forecasting without refitting. The observed
%   validation-year multivariate history only updates the test input.
%   No future exogenous observations enter a forecast.
%
% Evaluation
%   NRMSE = RMSE / (max(test dGWS) - min(test dGWS)).
%   TDA compares signs of successive monthly dGWS differences, using
%   11 within-year comparisons per test year.
%
% Output
%   RFX_Pred, RFX_Val_Pred: [lat, lon, 36].
%   RFX_NRMSE, RFX_TDA, RFX_best_p, RFX_best_x,
%   RFX_Val_NRMSE_best: [lat, lon, 3].
%   RFX_RMSE, RFX_TestRange, RFX_TestSuccess: [lat, lon, 3].
%   RFX_FailureCode: [lat, lon, 3]; code definitions in Meta.
%   Valid-grid arrays, spatial coordinates, dates and metadata are saved.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
% Pilot wrapper: do not clear the input configuration.
projectRoot = cfg.projectRoot;
dataFile = cfg.pilotMasterFile;
outFolder = fullfile(cfg.outputRoot, 'RFX');

% Independent experiments with distinct random seeds and output files.
experimentName = 'original';  % 'original', 'uc1', or 'uc2'
switch lower(experimentName)
    case 'original'
        randomSeed = 20260605;
        outFileName = 'RFX_extended_2019_2021.mat';
        checkpointFileName = 'RFX_extended_2019_2021_checkpoint.mat';
    case 'uc1'
        randomSeed = 20260719;
        outFileName = 'RFX_multivariate_global_200401_202412_uc1.mat';
        checkpointFileName = 'RFX_multivariate_global_uc1_checkpoint.mat';
    case 'uc2'
        randomSeed = 20260721;
        outFileName = 'RFX_multivariate_global_200401_202412_uc2.mat';
        checkpointFileName = 'RFX_multivariate_global_uc2_checkpoint.mat';
    otherwise
        error('experimentName must be ''original'', ''uc1'', or ''uc2''.');
end

% Existing final output: "overwrite" or "stop".
outputMode = "overwrite";
resumeFromCheckpoint = true;

useParallel = cfg.useParallel;
numWorkers = cfg.numWorkers;
block_num = cfg.block_num;

p_list = [1, 6, 12, 24, 36, 48, 60];
x_list = [1, 3, 6];
inputVarNames = {'dGWS', 'Prec', 'Evap', 'Runoff', 'dSM2m', 'SurfT'};
exogVarNames = {'Prec', 'Evap', 'Runoff', 'dSM2m', 'SurfT'};
testYears = [2019, 2020, 2021];
forecastHorizon = 12;

rfParam = struct();
rfParam.Method = 'Bag';
rfParam.NumLearningCycles = 50;
rfParam.MinLeafSize = 5;
rfParam.NumVariablesToSampleMode = 'p_x_over_3';
numRepeat = 1;

makeDiagnosticPlot = false;
plotNRMSELimits = [0, 0.4];
checkpointVersion = 2;

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
if any(x_list < 1) || any(x_list ~= floor(x_list)) || ...
        numel(unique(x_list)) ~= numel(x_list)
    error('x_list must contain distinct positive integers.');
end
if ~isscalar(randomSeed) || randomSeed < 1 || ...
        randomSeed ~= floor(randomSeed) || randomSeed > 2^32 - 1
    error('randomSeed must be a positive integer within the uint32 range.');
end
if ~isequal(testYears, [2019, 2020, 2021]) || forecastHorizon ~= 12
    error('This experiment requires three annual 12-month forecasts, 2019-2021.');
end
if ~isscalar(numRepeat) || numRepeat < 1 || numRepeat ~= floor(numRepeat)
    error('numRepeat must be a positive integer.');
end
if ~strcmpi(rfParam.Method, 'Bag') || ...
        ~isscalar(rfParam.NumLearningCycles) || ...
        rfParam.NumLearningCycles < 1 || ...
        rfParam.NumLearningCycles ~= floor(rfParam.NumLearningCycles) || ...
        ~isscalar(rfParam.MinLeafSize) || rfParam.MinLeafSize < 1 || ...
        rfParam.MinLeafSize ~= floor(rfParam.MinLeafSize) || ...
        ~ismember(rfParam.NumVariablesToSampleMode, {'p_x_over_3','sqrt','all'})
    error('Invalid RFX ensemble settings.');
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

requiredVars = {'dGWS_valid', 'Prec_valid', 'Evap_valid', ...
    'Runoff_valid', 'dSM2m_valid', 'SurfT_valid', ...
    'CommonMask', 'CommonValid', 'validLinearIndex', 'validLat', ...
    'validLon', 'lat', 'lon', 'ym'};
requiredVars{end+1} = 'originalValidRows';
checkVarsExist(who('-file', dataFile), requiredVars);

D = load(dataFile, requiredVars{:});
dGWS_valid = double(D.dGWS_valid);
Prec_valid = double(D.Prec_valid);
Evap_valid = double(D.Evap_valid);
Runoff_valid = double(D.Runoff_valid);
dSM2m_valid = double(D.dSM2m_valid);
SurfT_valid = double(D.SurfT_valid);
CommonMask = D.CommonMask;
CommonValid = logical(D.CommonValid);
validLinearIndex = double(D.validLinearIndex(:));
validLat = double(D.validLat(:));
validLon = double(D.validLon(:));
lat = double(D.lat(:));
lon = double(D.lon(:));
ym = double(D.ym(:));
originalValidRows = double(D.originalValidRows(:));
clear D

nLat = numel(lat);
nLon = numel(lon);
nValid = size(dGWS_valid, 1);
nTime = numel(ym);
expectedYM = makeYMList(2004, 1, 2024, 12);

if nLat ~= 150 || nLon ~= 360 || ~isequal(ym, expectedYM)
    error('Expected a 150 x 360 grid and months from 2004-01 to 2024-12.');
end
inputData = {dGWS_valid, Prec_valid, Evap_valid, Runoff_valid, ...
    dSM2m_valid, SurfT_valid};
for v = 1:numel(inputData)
    if ~isequal(size(inputData{v}), [nValid, nTime]) || ...
            any(~isfinite(inputData{v}(:)))
        error('%s must contain finite [grid, 252] time series.', ...
            inputVarNames{v});
    end
end
clear inputData
if nTime ~= 252
    error('Expected 252 monthly observations.');
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
    error('Valid-grid coordinates do not match validLinearIndex.');
end

if numel(originalValidRows) ~= nValid || ...
        any(~isfinite(originalValidRows)) || ...
        any(originalValidRows < 1) || ...
        any(originalValidRows ~= floor(originalValidRows)) || ...
        any(diff(originalValidRows) <= 0)
    error('Pilot source valid-grid row IDs are inconsistent.');
end

nRound = numel(testYears);
nPredMonth = nRound * forecastHorizon;
predYM = makeYMList(testYears(1), 1, testYears(end), 12);
valYM = makeYMList(testYears(1)-1, 1, testYears(end)-1, 12);
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);

block_size = max(1, ceil(nValid / block_num));
total_blocks = ceil(nValid / block_size);
fprintf('RFX experiment: %s, random seed %d.\n', experimentName, randomSeed);
fprintf('Grid cells: %d; monthly observations: %d; test rounds: %d.\n', ...
    nValid, nTime, nRound);
fprintf('History lengths: %s\n', mat2str(p_list));
fprintf('Exogenous history lengths: %s\n', mat2str(x_list));
fprintf('Bagged trees: %d; leaf size: %d; predictor sampling: %s.\n', ...
    rfParam.NumLearningCycles, rfParam.MinLeafSize, ...
    rfParam.NumVariablesToSampleMode);
fprintf('Replicates per candidate: %d.\n', numRepeat);
fprintf('Blocks: %d; approximately %d grid cells per block.\n', ...
    total_blocks, block_size);

%% 3.2 Prepare checkpoint and result arrays
runSignature = computeRunSignature(dGWS_valid, Prec_valid, Evap_valid, ...
    Runoff_valid, dSM2m_valid, SurfT_valid, validLinearIndex, ...
    lat, lon, ym, p_list, x_list, testYears, forecastHorizon, ...
    numRepeat, rfParam, randomSeed, block_num, checkpointVersion);

if isfile(checkpointFile) && ~resumeFromCheckpoint
    delete(checkpointFile);
end

if resumeFromCheckpoint && isfile(checkpointFile)
    fprintf('Reading checkpoint: %s\n', checkpointFile);
    C = load(checkpointFile);
    validateCheckpoint(C, runSignature, nValid, nRound, nPredMonth, total_blocks);

    RFX_Pred_valid = C.RFX_Pred_valid;
    RFX_Val_Pred_valid = C.RFX_Val_Pred_valid;
    RFX_NRMSE_valid = C.RFX_NRMSE_valid;
    RFX_TDA_valid = C.RFX_TDA_valid;
    RFX_best_p_valid = C.RFX_best_p_valid;
    RFX_best_x_valid = C.RFX_best_x_valid;
    RFX_Val_NRMSE_best_valid = C.RFX_Val_NRMSE_best_valid;
    RFX_RMSE_valid = C.RFX_RMSE_valid;
    RFX_TestRange_valid = C.RFX_TestRange_valid;
    RFX_TestSuccess_valid = C.RFX_TestSuccess_valid;
    RFX_FailureCode_valid = C.RFX_FailureCode_valid;
    RFX_FailureID_valid = C.RFX_FailureID_valid;
    blockDone = C.blockDone;
    cumulativeComputeTimeSeconds = C.cumulativeComputeTimeSeconds;
    clear C
else
    RFX_Pred_valid = NaN(nValid, nPredMonth);
    RFX_Val_Pred_valid = NaN(nValid, nPredMonth);
    RFX_NRMSE_valid = NaN(nValid, nRound);
    RFX_TDA_valid = NaN(nValid, nRound);
    RFX_best_p_valid = NaN(nValid, nRound);
    RFX_best_x_valid = NaN(nValid, nRound);
    RFX_Val_NRMSE_best_valid = NaN(nValid, nRound);
    RFX_RMSE_valid = NaN(nValid, nRound);
    RFX_TestRange_valid = NaN(nValid, nRound);
    RFX_TestSuccess_valid = false(nValid, nRound);
    RFX_FailureCode_valid = uint8(zeros(nValid, nRound));
    RFX_FailureID_valid = cell(nValid, nRound);
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
    Prec_const = parallel.pool.Constant(Prec_valid);
    Evap_const = parallel.pool.Constant(Evap_valid);
    Runoff_const = parallel.pool.Constant(Runoff_valid);
    dSM2m_const = parallel.pool.Constant(dSM2m_valid);
    SurfT_const = parallel.pool.Constant(SurfT_valid);
end

%% 3.4 Run global RFX forecasts
fprintf('\nStarting RFX computation...\n');
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
    blockBestX = NaN(block_len, nRound);
    blockValBest = NaN(block_len, nRound);
    blockRMSE = NaN(block_len, nRound);
    blockTestRange = NaN(block_len, nRound);
    blockSuccess = false(block_len, nRound);
    blockFailureCode = uint8(zeros(block_len, nRound));
    blockFailureID = cell(block_len, nRound);

    blockTic = tic;
    if useParallel
        parfor i = 1:block_len
            gridRow = currentRows(i);
            X_pixel = [Y_valid_const.Value(gridRow, :)', ...
                Prec_const.Value(gridRow, :)', ...
                Evap_const.Value(gridRow, :)', ...
                Runoff_const.Value(gridRow, :)', ...
                dSM2m_const.Value(gridRow, :)', ...
                SurfT_const.Value(gridRow, :)'];
            [p36, v36, n3, t3, best3, x3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runRFXOnePixelAllRounds(X_pixel, ...
                p_list, x_list, testYears, forecastHorizon, ...
                numRepeat, rfParam, randomSeed, originalValidRows(gridRow));
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockBestX(i, :) = x3(:)';
            blockValBest(i, :) = vn3(:)';
            blockRMSE(i, :) = rmse3(:)';
            blockTestRange(i, :) = range3(:)';
            blockSuccess(i, :) = ok3(:)';
            blockFailureCode(i, :) = code3(:)';
            blockFailureID(i, :) = id3(:)';
        end
    else
        for i = 1:block_len
            gridRow = currentRows(i);
            X_pixel = [dGWS_valid(gridRow, :)', ...
                Prec_valid(gridRow, :)', ...
                Evap_valid(gridRow, :)', ...
                Runoff_valid(gridRow, :)', ...
                dSM2m_valid(gridRow, :)', ...
                SurfT_valid(gridRow, :)'];
            [p36, v36, n3, t3, best3, x3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runRFXOnePixelAllRounds(X_pixel, ...
                p_list, x_list, testYears, forecastHorizon, ...
                numRepeat, rfParam, randomSeed, originalValidRows(gridRow));
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockBestX(i, :) = x3(:)';
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

    RFX_Pred_valid(currentRows, :) = blockPred;
    RFX_Val_Pred_valid(currentRows, :) = blockValPred;
    RFX_NRMSE_valid(currentRows, :) = blockNRMSE;
    RFX_TDA_valid(currentRows, :) = blockTDA;
    RFX_best_p_valid(currentRows, :) = blockBestP;
    RFX_best_x_valid(currentRows, :) = blockBestX;
    RFX_Val_NRMSE_best_valid(currentRows, :) = blockValBest;
    RFX_RMSE_valid(currentRows, :) = blockRMSE;
    RFX_TestRange_valid(currentRows, :) = blockTestRange;
    RFX_TestSuccess_valid(currentRows, :) = blockSuccess;
    RFX_FailureCode_valid(currentRows, :) = blockFailureCode;
    RFX_FailureID_valid(currentRows, :) = blockFailureID;
    blockDone(block) = true;

    checkpointInfo = struct('version', checkpointVersion, ...
        'runSignature', runSignature);
    checkpointTemp = fullfile(outFolder, [checkpointFileName, '.tmp.mat']);
    save(checkpointTemp, ...
        'RFX_Pred_valid', 'RFX_Val_Pred_valid', ...
        'RFX_NRMSE_valid', 'RFX_TDA_valid', ...
        'RFX_best_p_valid', 'RFX_best_x_valid', ...
        'RFX_Val_NRMSE_best_valid', 'RFX_RMSE_valid', 'RFX_TestRange_valid', ...
        'RFX_TestSuccess_valid', 'RFX_FailureCode_valid', 'RFX_FailureID_valid', ...
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
fprintf('\nRFX computation finished.\n');
fprintf('Computation time this run: %.2f minutes (%.2f hours)\n', ...
    elapsedComputeThisRun / 60, elapsedComputeThisRun / 3600);
fprintf('Cumulative block computation time: %.2f minutes (%.2f hours)\n', ...
    cumulativeComputeTimeSeconds / 60, cumulativeComputeTimeSeconds / 3600);
if useParallel
    delete(Y_valid_const);
    delete(Prec_const);
    delete(Evap_const);
    delete(Runoff_const);
    delete(dSM2m_const);
    delete(SurfT_const);
end

%% 3.5 Map valid-grid results to 1-degree fields
RFX_Pred = mapValidToGrid3D(RFX_Pred_valid, validLinearIndex, nLat, nLon);
RFX_Val_Pred = mapValidToGrid3D(RFX_Val_Pred_valid, validLinearIndex, nLat, nLon);
RFX_NRMSE = mapValidToGrid3D(RFX_NRMSE_valid, validLinearIndex, nLat, nLon);
RFX_TDA = mapValidToGrid3D(RFX_TDA_valid, validLinearIndex, nLat, nLon);
RFX_best_p = mapValidToGrid3D(RFX_best_p_valid, validLinearIndex, nLat, nLon);
RFX_best_x = mapValidToGrid3D(RFX_best_x_valid, validLinearIndex, nLat, nLon);
RFX_Val_NRMSE_best = mapValidToGrid3D(RFX_Val_NRMSE_best_valid, ...
    validLinearIndex, nLat, nLon);
RFX_RMSE = mapValidToGrid3D(RFX_RMSE_valid, validLinearIndex, nLat, nLon);
RFX_TestRange = mapValidToGrid3D(RFX_TestRange_valid, validLinearIndex, nLat, nLon);
RFX_TestSuccess = mapValidToGrid3D(double(RFX_TestSuccess_valid), ...
    validLinearIndex, nLat, nLon);
RFX_FailureCode = mapValidToGrid3D(double(RFX_FailureCode_valid), ...
    validLinearIndex, nLat, nLon);

RFX_NRMSE_mean = mean(RFX_NRMSE, 3, 'omitnan');
RFX_TDA_mean = mean(RFX_TDA, 3, 'omitnan');
RFX_Val_NRMSE_best_mean = mean(RFX_Val_NRMSE_best, 3, 'omitnan');

CompleteThreeYears_valid = all(RFX_TestSuccess_valid, 2);
RFX_CompleteThreeYears = false(nLat, nLon);
RFX_CompleteThreeYears(validLinearIndex) = CompleteThreeYears_valid;
RFX_NRMSE_mean_complete = RFX_NRMSE_mean;
RFX_NRMSE_mean_complete(~RFX_CompleteThreeYears) = NaN;

fprintf('\nRFX global results:\n');
fprintf('Successful grid-year forecasts: %d / %d.\n', ...
    nnz(RFX_TestSuccess_valid), numel(RFX_TestSuccess_valid));
fprintf('Complete three-year grid cells: %d / %d.\n', ...
    nnz(CompleteThreeYears_valid), nValid);
for r = 1:nRound
    fprintf('  %d success: %d / %d\n', testYears(r), ...
        nnz(RFX_TestSuccess_valid(:, r)), nValid);
    for code = 1:4
        nFailure = nnz(RFX_FailureCode_valid(:, r) == code);
        if nFailure > 0
            fprintf('    Failure code %d: %d\n', code, nFailure);
        end
    end
end
fprintf('RFX_NRMSE_mean (finite values): %.4f\n', ...
    mean(RFX_NRMSE_mean(CommonValid), 'omitnan'));
fprintf('RFX_TDA_mean (finite values): %.2f%%\n', ...
    mean(RFX_TDA_mean(CommonValid), 'omitnan'));

%% 3.6 Save final output
Meta = struct();
Meta.model = 'RFX';
Meta.experimentName = experimentName;
Meta.input = ['Historical dGWS plus precipitation, evaporation, runoff, ' ...
    'soil-moisture change and surface temperature'];
Meta.inputVarNames = inputVarNames;
Meta.exogVarNames = exogVarNames;
Meta.target = 'Monthly dGWS';
Meta.period = '2004-01 to 2024-12';
Meta.predictionPeriod = '2019-01 to 2021-12';
Meta.validationPeriod = '2018-01 to 2020-12';
Meta.arrayConvention = 'Maps: [latitude, longitude, time]; valid series: [grid, time]';
Meta.rollingStrategy = ['Train 2004 through test year minus two; ' ...
    'validate in the following year; test in the subsequent year'];
Meta.p_list = p_list;
Meta.x_list = x_list;
Meta.numRepeat = numRepeat;
Meta.rfParam = rfParam;
Meta.randomSeed = randomSeed;
Meta.standardization = 'None; bagged regression trees use original lagged values';
Meta.forecasting = 'Twelve direct bagged-tree ensembles, one for each lead month';
Meta.noFutureExogenousVariables = ['Only lagged exogenous variables that ' ...
    'are observed before each forecast origin are used'];
Meta.finalPrediction = ['Validation-selected fitted RF ensembles are reused ' ...
    'for testing without refitting; observed validation-year multivariate ' ...
    'history updates only predictors at the test forecast origin'];
Meta.repeatStrategy = ['The validation score is the mean of finite replicate-specific ' ...
    'NRMSE values. Validation and test predictions average selected fitted ' ...
    'ensembles across successful replicates; the averaged validation NRMSE ' ...
    'is not necessarily the selection score when numRepeat exceeds one'];
Meta.randomization = ['Each ensemble fit uses a seed based on experiment, ' ...
    'grid row, test round, p, x, replicate and horizon. ' ...
    'Predictor selection is reproducible and no test-stage fit is made'];
Meta.NRMSE = 'Annual RMSE divided by observed test-year dGWS range';
Meta.TDA = ['Percentage of matching signs of successive month-to-month ' ...
    'differences in dGWS within a test year (11 comparisons)'];
Meta.failureCodes = ['0 = success; 1 = no finite validation candidate; ' ...
    '2 = test forecast error; 3 = invalid test predictions; ' ...
    '4 = nonfinite test NRMSE'];
Meta.RFX_FailureID = 'MATLAB exception identifiers for failed fits or forecasts';
Meta.meanMetricRule = 'Three-year mean computed with omitnan';
Meta.completeMetricRule = ['RFX_NRMSE_mean_complete requires successful ' ...
    'forecasts in all three test years'];
Meta.nValid = nValid;
Meta.block_num = block_num;
Meta.numWorkers = numWorkers;
Meta.usedParallel = useParallel;
Meta.runSignature = runSignature;
Meta.computeTimeMinutes = cumulativeComputeTimeSeconds / 60;
Meta.computeTimeThisRunMinutes = elapsedComputeThisRun / 60;
Meta.scriptName = 'SKLGB_ext_RFX.m';
Meta.pilotTag = 'extended_2019_2021';
Meta.pilotOriginalValidRows = originalValidRows;
Meta.pilotN = nValid;
Meta.pilotMode = 'extended-2019-2021-original-hyperparameters';

finalTemp = fullfile(outFolder, [outFileName, '.tmp.mat']);
fprintf('\nSaving final MAT-file: %s\n', outFile);
save(finalTemp, ...
    'RFX_Pred', 'RFX_Val_Pred', ...
    'RFX_NRMSE', 'RFX_NRMSE_mean', 'RFX_NRMSE_mean_complete', ...
    'RFX_TDA', 'RFX_TDA_mean', 'RFX_best_p', 'RFX_best_x', ...
    'RFX_Val_NRMSE_best', 'RFX_Val_NRMSE_best_mean', ...
    'RFX_RMSE', 'RFX_TestRange', 'RFX_TestSuccess', ...
    'RFX_FailureCode', 'RFX_CompleteThreeYears', ...
    'RFX_Pred_valid', 'RFX_Val_Pred_valid', ...
    'RFX_NRMSE_valid', 'RFX_TDA_valid', 'RFX_best_p_valid', ...
    'RFX_best_x_valid', ...
    'RFX_Val_NRMSE_best_valid', ...
    'RFX_RMSE_valid', 'RFX_TestRange_valid', 'RFX_TestSuccess_valid', ...
    'RFX_FailureCode_valid', 'RFX_FailureID_valid', ...
    'p_list', 'x_list', 'numRepeat', 'testYears', 'predYM', 'valYM', ...
    'rfParam', 'randomSeed', 'experimentName', ...
    'inputVarNames', 'exogVarNames', ...
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
    figure('Name', 'RFX mean NRMSE');
    p = imagesc(lon, lat, RFX_NRMSE_mean);
    set(p, 'AlphaData', double(isfinite(RFX_NRMSE_mean)));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(jet);
    clim(plotNRMSELimits);
    cb = colorbar;
    ylabel(cb, 'NRMSE');
    title('RFX mean NRMSE, test years 2019-2021');
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

elapsedTotalSeconds = toc(totalTic);
fprintf('\nRFX processing completed.\n');
fprintf('Total time since computation started (including output): %.2f minutes (%.2f hours)\n', ...
    elapsedTotalSeconds / 60, elapsedTotalSeconds / 3600);

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

end  % End of SKLGB_pilot_RFX main function

function [pred36, valPred36, nrmse3, tda3, bestP3, bestX3, ...
    valBest3, rmse3, testRange3, success3, failureCode3, failureID3] = ...
    runRFXOnePixelAllRounds(X_all, p_list, x_list, testYears, horizon, ...
    numRepeat, rfParam, randomSeed, pointID)

    X_all = double(X_all);
    Y = X_all(:, 1);
    nRound = numel(testYears);
    pred36 = NaN(nRound * horizon, 1);
    valPred36 = NaN(nRound * horizon, 1);
    nrmse3 = NaN(nRound, 1);
    tda3 = NaN(nRound, 1);
    bestP3 = NaN(nRound, 1);
    bestX3 = NaN(nRound, 1);
    valBest3 = NaN(nRound, 1);
    rmse3 = NaN(nRound, 1);
    testRange3 = NaN(nRound, 1);
    success3 = false(nRound, 1);
    failureCode3 = uint8(zeros(nRound, 1));
    failureID3 = cell(nRound, 1);

    if ~isequal(size(X_all), [252, 6]) || any(~isfinite(X_all(:)))
        error('RFX input must contain six finite, 252-month time series.');
    end

    for r = 1:nRound
        testYear = testYears(r);
        trainEndIdx = (testYear - 2005) * 12;
        valIdx = (trainEndIdx + 1):(trainEndIdx + horizon);
        testIdx = (trainEndIdx + horizon + 1):(trainEndIdx + 2 * horizon);
        predIdx = ((r - 1) * horizon + 1):(r * horizon);

        X_train = X_all(1:trainEndIdx, :);
        Y_val = Y(valIdx);
        Y_test = Y(testIdx);
        bestValScore = Inf;
        bestP = NaN;
        bestX = NaN;
        selectedPacks = cell(numRepeat, 1);
        selectedValPred = NaN(horizon, 1);
        firstValidationError = '';

        for ip = 1:numel(p_list)
            p = p_list(ip);
            for ix = 1:numel(x_list)
                x = x_list(ix);
                maxLag = max(p, x);
                if size(X_train, 1) <= maxLag + horizon - 1
                    continue;
                end

                valScores = NaN(numRepeat, 1);
                candidatePacks = cell(numRepeat, 1);
                candidateValPred = NaN(horizon, numRepeat);
                initWindowVal = X_train(end-maxLag+1:end, :);
                for rep = 1:numRepeat
                    try
                        seedVal = makeSeed(randomSeed, pointID, r, ip, ix, rep);
                        pack = fitRFXHorizonSpecificFromSeries( ...
                            X_train, p, x, horizon, rfParam, seedVal);
                        predVal = predictRFXHorizonSpecific(pack, initWindowVal);
                        score = calcNRMSE(Y_val, predVal);
                        if isfinite(score) && numel(predVal) == horizon && ...
                                all(isfinite(predVal))
                            valScores(rep) = score;
                            candidatePacks{rep} = pack;
                            candidateValPred(:, rep) = predVal(:);
                        end
                    catch ME
                        if isempty(firstValidationError)
                            firstValidationError = ME.identifier;
                        end
                    end
                end

                validRep = isfinite(valScores);
                if ~any(validRep)
                    continue;
                end
                candidateScore = mean(valScores(validRep));
                if isfinite(candidateScore) && candidateScore < bestValScore
                    bestValScore = candidateScore;
                    bestP = p;
                    bestX = x;
                    selectedPacks = candidatePacks;
                    selectedValPred = mean(candidateValPred(:, validRep), 2);
                end
            end
        end

        if ~isfinite(bestValScore)
            failureCode3(r) = 1;
            failureID3{r} = firstValidationError;
            continue;
        end

        bestP3(r) = bestP;
        bestX3(r) = bestX;
        valBest3(r) = bestValScore;
        valPred36(predIdx) = selectedValPred(:);

        % Only the forecast-origin history moves forward by one year.
        maxLag = max(bestP, bestX);
        originIdx = trainEndIdx + horizon;
        initWindowTest = X_all(originIdx-maxLag+1:originIdx, :);
        testPredRepeats = NaN(horizon, numRepeat);
        firstTestError = '';
        for rep = 1:numRepeat
            pack = selectedPacks{rep};
            if isempty(pack)
                continue;
            end
            try
                testPredRepeats(:, rep) = predictRFXHorizonSpecific( ...
                    pack, initWindowTest);
            catch ME
                if isempty(firstTestError)
                    firstTestError = ME.identifier;
                end
            end
        end

        okRep = all(isfinite(testPredRepeats), 1);
        if ~any(okRep)
            failureCode3(r) = 2;
            failureID3{r} = firstTestError;
            continue;
        end
        Y_test_pred = mean(testPredRepeats(:, okRep), 2);
        if numel(Y_test_pred) ~= horizon || any(~isfinite(Y_test_pred))
            failureCode3(r) = 3;
            continue;
        end
        pred36(predIdx) = Y_test_pred(:);
        testNRMSE = calcNRMSE(Y_test, Y_test_pred);
        if ~isfinite(testNRMSE)
            failureCode3(r) = 4;
            continue;
        end
        nrmse3(r) = testNRMSE;
        rmse3(r) = sqrt(mean((Y_test_pred(:) - Y_test).^2));
        testRange3(r) = max(Y_test) - min(Y_test);
        tda3(r) = calcTDA(Y_test, Y_test_pred);
        success3(r) = true;
    end
end

function rfPack = fitRFXHorizonSpecificFromSeries( ...
    X_all, p, x, horizon, rfParam, baseSeed)

    X_all = double(X_all);
    nPredictor = p + (size(X_all, 2) - 1) * x;
    numVariablesToSample = getNumVariablesToSample(nPredictor, rfParam);
    treeTemplate = templateTree( ...
        'MinLeafSize', rfParam.MinLeafSize, ...
        'NumVariablesToSample', numVariablesToSample, ...
        'Reproducible', true);
    modelCell = cell(horizon, 1);
    for h = 1:horizon
        [X_train, Y_target] = buildRFXLaggedSamplesDirect(X_all, p, x, h);
        seedH = mod(baseSeed + h, 2^32 - 1);
        if seedH == 0
            seedH = 1;
        end
        rng(seedH, 'twister');
        modelCell{h} = fitrensemble(X_train, Y_target, ...
            'Method', rfParam.Method, ...
            'NumLearningCycles', rfParam.NumLearningCycles, ...
            'Learners', treeTemplate);
    end
    rfPack = struct();
    rfPack.modelCell = modelCell;
    rfPack.p = p;
    rfPack.x = x;
    rfPack.horizon = horizon;
    rfPack.nInputVars = size(X_all, 2);
end

function [X_train, Y_target] = buildRFXLaggedSamplesDirect(X_all, p, x, h)
    X_all = double(X_all);
    nTime = size(X_all, 1);
    nInputVars = size(X_all, 2);
    nExogVars = nInputVars - 1;
    maxLag = max(p, x);
    nSample = nTime - maxLag - h + 1;
    if nSample <= 0
        error('No RFX direct samples can be constructed.');
    end
    X_train = zeros(nSample, p + nExogVars * x);
    Y_target = zeros(nSample, 1);
    for j = 1:nSample
        anchorIdx = j + maxLag - 1;
        dGWSLag = X_all(anchorIdx-p+1:anchorIdx, 1)';
        exogFeat = zeros(1, nExogVars * x);
        for v = 2:nInputVars
            featurePos = (v-2)*x + (1:x);
            exogFeat(featurePos) = X_all(anchorIdx-x+1:anchorIdx, v)';
        end
        X_train(j, :) = [dGWSLag, exogFeat];
        Y_target(j) = X_all(anchorIdx+h, 1);
    end
end

function Y_pred = predictRFXHorizonSpecific(rfPack, initialWindow)
    p = rfPack.p;
    x = rfPack.x;
    horizon = rfPack.horizon;
    maxLag = max(p, x);
    initialWindow = double(initialWindow);
    if ~isequal(size(initialWindow), [maxLag, rfPack.nInputVars])
        error('The RFX initial history does not match p, x and variables.');
    end
    predictor = makeRFXOneFeatureRow(initialWindow, p, x);
    Y_pred = NaN(horizon, 1);
    for h = 1:horizon
        prediction = predict(rfPack.modelCell{h}, predictor);
        Y_pred(h) = double(prediction(1));
    end
end

function row = makeRFXOneFeatureRow(initialWindow, p, x)
    nInputVars = size(initialWindow, 2);
    nExogVars = nInputVars - 1;
    maxLag = max(p, x);
    anchorIdx = maxLag;
    row = zeros(1, p + nExogVars*x);
    row(1:p) = initialWindow(anchorIdx-p+1:anchorIdx, 1)';
    for v = 2:nInputVars
        pos = p + (v-2)*x + (1:x);
        row(pos) = initialWindow(anchorIdx-x+1:anchorIdx, v)';
    end
end

function nVars = getNumVariablesToSample(nPredictor, rfParam)
    switch rfParam.NumVariablesToSampleMode
        case 'p_x_over_3'
            nVars = max(1, floor(nPredictor / 3));
        case 'sqrt'
            nVars = max(1, floor(sqrt(nPredictor)));
        case 'all'
            nVars = nPredictor;
        otherwise
            error('Unknown NumVariablesToSampleMode.');
    end
    nVars = min(nVars, nPredictor);
end

function seedVal = makeSeed(randomSeed, pointID, roundID, pID, xID, repeatID)
    seedRaw = double(randomSeed) + double(pointID)*100 + ...
        double(roundID)*10000 + double(pID)*1000 + ...
        double(xID)*100 + double(repeatID)*10;
    seedVal = mod(seedRaw, 2^32 - 1);
    if seedVal == 0
        seedVal = randomSeed;
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
    tda = mean(sign(diff(obs)) == sign(diff(pred))) * 100;
end

function Ygrid = mapValidToGrid3D(Yvalid, validLinearIndex, nLat, nLon)
    [~, nLayer] = size(Yvalid);
    Ygrid = NaN(nLat, nLon, nLayer);
    for k = 1:nLayer
        layer = NaN(nLat, nLon);
        layer(validLinearIndex) = Yvalid(:, k);
        Ygrid(:, :, k) = layer;
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

function signature = computeRunSignature(dGWS, prec, evap, runoff, dSM2m, ...
    surfT, index, lat, lon, ym, p_list, x_list, testYears, horizon, ...
    numRepeat, rfParam, randomSeed, blockNum, version)
    md = java.security.MessageDigest.getInstance('SHA-256');
    updateDigest(md, 'RFX_global_forecast');
    updateDigest(md, double(version));
    updateDigest(md, dGWS);
    updateDigest(md, prec);
    updateDigest(md, evap);
    updateDigest(md, runoff);
    updateDigest(md, dSM2m);
    updateDigest(md, surfT);
    updateDigest(md, index);
    updateDigest(md, lat);
    updateDigest(md, lon);
    updateDigest(md, ym);
    updateDigest(md, double(p_list));
    updateDigest(md, double(x_list));
    updateDigest(md, double(testYears));
    updateDigest(md, double(horizon));
    updateDigest(md, double(numRepeat));
    updateDigest(md, rfParam.Method);
    updateDigest(md, double(rfParam.NumLearningCycles));
    updateDigest(md, double(rfParam.MinLeafSize));
    updateDigest(md, rfParam.NumVariablesToSampleMode);
    updateDigest(md, double(randomSeed));
    updateDigest(md, double(blockNum));
    hashBytes = typecast(int8(md.digest()), 'uint8');
    signature = lower(reshape(dec2hex(hashBytes, 2).', 1, []));
end

function updateDigest(md, value)
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
        'RFX_Pred_valid', 'RFX_Val_Pred_valid', ...
        'RFX_NRMSE_valid', 'RFX_TDA_valid', ...
        'RFX_best_p_valid', 'RFX_best_x_valid', ...
        'RFX_Val_NRMSE_best_valid', 'RFX_RMSE_valid', ...
        'RFX_TestRange_valid', 'RFX_TestSuccess_valid', ...
        'RFX_FailureCode_valid', 'RFX_FailureID_valid', ...
        'cumulativeComputeTimeSeconds'};
    for k = 1:numel(required)
        if ~isfield(C, required{k})
            error('Checkpoint is missing %s. Use resumeFromCheckpoint=false.', ...
                required{k});
        end
    end
    if ~isfield(C.checkpointInfo, 'runSignature') || ...
            ~strcmp(C.checkpointInfo.runSignature, signature)
        error(['Checkpoint input data or settings differ from this run. ' ...
            'Use resumeFromCheckpoint=false to restart.']);
    end
    matrix3 = {'RFX_NRMSE_valid', 'RFX_TDA_valid', ...
        'RFX_best_p_valid', 'RFX_best_x_valid', ...
        'RFX_Val_NRMSE_best_valid', 'RFX_RMSE_valid', ...
        'RFX_TestRange_valid', 'RFX_TestSuccess_valid', ...
        'RFX_FailureCode_valid', 'RFX_FailureID_valid'};
    for k = 1:numel(matrix3)
        if ~isequal(size(C.(matrix3{k})), [nValid, nRound])
            error('Checkpoint contains an invalid %s array.', matrix3{k});
        end
    end
    if ~isequal(size(C.RFX_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.RFX_Val_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.blockDone), [totalBlocks, 1])
        error('Checkpoint dimensions do not match the current run.');
    end
end
