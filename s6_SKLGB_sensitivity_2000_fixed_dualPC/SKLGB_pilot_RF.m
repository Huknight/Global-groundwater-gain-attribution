function SKLGB_pilot_RF(cfg)
% PILOT WRAPPER: adapted from RF_global_forecast.m.
% Only file paths, pilot seed IDs, run settings and the two tested
% hyperparameters are changed; the original model code follows.
% Each configuration searches the FULL original validation p/x grid.

%% RF - Global univariate groundwater forecasting
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Forecast monthly groundwater-storage change (dGWS) for each valid
%   1-degree grid cell using bagged regression trees.
%
% Input
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   dGWS_valid: [grid, month], January 2004 through December 2024.
%
% Rolling evaluation
%   Test 2022: train 2004-2020, validate 2021.
%   Test 2023: train 2004-2021, validate 2022.
%   Test 2024: train 2004-2022, validate 2023.
%   Each test-year forecast is issued from the preceding December.
%
% Model
%   Bagged regression-tree ensemble, fitrensemble, 50 trees, leaf size 5.
%   History window p: [1, 6, 12, 24, 36, 48, 60] months.
%   Predictors sampled at each split: max(1, floor(p / 3)).
%   For each p, twelve separate ensembles forecast h = 1,...,12.
%   Validation-year NRMSE selects the history window length.
%   For numRepeat > 1, the selection score is the mean of the valid
%   replicate-specific validation NRMSE values; test predictions average
%   the corresponding selected fitted ensembles across replicates.
%   Validation-selected fitted ensembles are reused without refitting.
%   Observed validation-year values update the test forecast-origin history.
%
% Evaluation
%   NRMSE = RMSE / (max(test dGWS) - min(test dGWS)).
%   TDA is agreement in successive month-to-month dGWS change directions,
%   based on 11 within-year differences per test year.
%
% Output
%   RF_Pred, RF_Val_Pred: [lat, lon, 36], prediction months in order.
%   RF_NRMSE, RF_TDA: [lat, lon, 3], test years 2022-2024.
%   RF_best_p, RF_Val_NRMSE_best: [lat, lon, 3].
%   RF_RMSE, RF_TestRange, RF_TestSuccess: [lat, lon, 3].
%   RF_FailureCode: [lat, lon, 3]; code definitions in Meta.
%   Additional valid-grid outputs, coordinates and metadata are saved.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
% Pilot wrapper: do not clear the input configuration.
projectRoot = cfg.projectRoot;
dataFile = cfg.pilotMasterFile;
outFolder = fullfile(cfg.outputRoot, 'RF', cfg.tag);

experimentName = 'original';  % 'original', 'uc1', or 'uc2'
switch lower(experimentName)
    case 'original'
        randomSeed = 20260605;
        outFileName = 'RF_univariate_global_200401_202412.mat';
        checkpointFileName = 'RF_univariate_global_checkpoint.mat';
    case 'uc1'
        randomSeed = 20260719;
        outFileName = 'RF_univariate_global_200401_202412_uc1.mat';
        checkpointFileName = 'RF_univariate_global_uc1_checkpoint.mat';
    case 'uc2'
        randomSeed = 20260721;
        outFileName = 'RF_univariate_global_200401_202412_uc2.mat';
        checkpointFileName = 'RF_univariate_global_uc2_checkpoint.mat';
    otherwise
        error('experimentName must be ''original'', ''uc1'', or ''uc2''.');
end

outputMode = "overwrite";  % 'overwrite' or 'stop'
resumeFromCheckpoint = true;

useParallel = cfg.useParallel;
numWorkers = cfg.numWorkers;
block_num = cfg.block_num;

p_list = [1, 6, 12, 24, 36, 48, 60];
testYears = [2022, 2023, 2024];
forecastHorizon = 12;
numRepeat = 1;

rfParam = struct();
rfParam.Method = 'Bag';
rfParam.NumLearningCycles = cfg.trees;
rfParam.MinLeafSize = cfg.leaf;
rfParam.NumVariablesToSampleMode = 'p_over_3';

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
if ~isscalar(numRepeat) || numRepeat < 1 || ...
        numRepeat ~= floor(numRepeat)
    error('numRepeat must be a positive integer.');
end
if ~isscalar(randomSeed) || randomSeed < 1 || ...
        randomSeed ~= floor(randomSeed) || randomSeed > 2^32 - 1
    error('randomSeed must be a positive integer within the uint32 range.');
end
if ~isequal(testYears, [2022, 2023, 2024]) || forecastHorizon ~= 12
    error('This experiment requires three annual 12-month forecasts, 2022-2024.');
end
if ~strcmpi(rfParam.Method, 'Bag') || ...
        ~isscalar(rfParam.NumLearningCycles) || ...
        rfParam.NumLearningCycles < 1 || ...
        rfParam.NumLearningCycles ~= floor(rfParam.NumLearningCycles) || ...
        ~isscalar(rfParam.MinLeafSize) || rfParam.MinLeafSize < 1 || ...
        rfParam.MinLeafSize ~= floor(rfParam.MinLeafSize) || ...
        ~ismember(rfParam.NumVariablesToSampleMode, {'p_over_3','sqrt','all'})
    error('Invalid RF ensemble settings.');
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
requiredVars{end+1} = 'originalValidRows';
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
fprintf('RF experiment: %s, random seed %d.\n', experimentName, randomSeed);
fprintf('Grid cells: %d; monthly observations: %d; test rounds: %d.\n', ...
    nValid, nTime, nRound);
fprintf('History lengths: %s\n', mat2str(p_list));
fprintf('Bagged trees: %d; leaf size: %d; predictor sampling: %s.\n', ...
    rfParam.NumLearningCycles, rfParam.MinLeafSize, ...
    rfParam.NumVariablesToSampleMode);
fprintf('Replicates per candidate: %d.\n', numRepeat);
fprintf('Blocks: %d; approximately %d grid cells per block.\n', ...
    total_blocks, block_size);

%% 3.2 Prepare checkpoint and result arrays
runSignature = computeRunSignature(dGWS_valid, validLinearIndex, ...
    lat, lon, ym, p_list, testYears, forecastHorizon, ...
    numRepeat, rfParam, randomSeed, block_num, checkpointVersion);

if isfile(checkpointFile) && ~resumeFromCheckpoint
    delete(checkpointFile);
end

if resumeFromCheckpoint && isfile(checkpointFile)
    fprintf('Reading checkpoint: %s\n', checkpointFile);
    C = load(checkpointFile);
    validateCheckpoint(C, runSignature, nValid, nRound, nPredMonth, total_blocks);

    RF_Pred_valid = C.RF_Pred_valid;
    RF_Val_Pred_valid = C.RF_Val_Pred_valid;
    RF_NRMSE_valid = C.RF_NRMSE_valid;
    RF_TDA_valid = C.RF_TDA_valid;
    RF_best_p_valid = C.RF_best_p_valid;
    RF_Val_NRMSE_best_valid = C.RF_Val_NRMSE_best_valid;
    RF_RMSE_valid = C.RF_RMSE_valid;
    RF_TestRange_valid = C.RF_TestRange_valid;
    RF_TestSuccess_valid = C.RF_TestSuccess_valid;
    RF_FailureCode_valid = C.RF_FailureCode_valid;
    RF_FailureID_valid = C.RF_FailureID_valid;
    blockDone = C.blockDone;
    cumulativeComputeTimeSeconds = C.cumulativeComputeTimeSeconds;
    clear C
else
    RF_Pred_valid = NaN(nValid, nPredMonth);
    RF_Val_Pred_valid = NaN(nValid, nPredMonth);
    RF_NRMSE_valid = NaN(nValid, nRound);
    RF_TDA_valid = NaN(nValid, nRound);
    RF_best_p_valid = NaN(nValid, nRound);
    RF_Val_NRMSE_best_valid = NaN(nValid, nRound);
    RF_RMSE_valid = NaN(nValid, nRound);
    RF_TestRange_valid = NaN(nValid, nRound);
    RF_TestSuccess_valid = false(nValid, nRound);
    RF_FailureCode_valid = uint8(zeros(nValid, nRound));
    RF_FailureID_valid = cell(nValid, nRound);
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

%% 3.4 Run global RF forecasts
fprintf('\nStarting RF computation...\n');
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
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runRFOnePixelAllRounds(Y_pixel, ...
                p_list, numRepeat, testYears, forecastHorizon, rfParam, ...
                randomSeed, originalValidRows(currentRows(i)));
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
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runRFOnePixelAllRounds(Y_pixel, ...
                p_list, numRepeat, testYears, forecastHorizon, rfParam, ...
                randomSeed, originalValidRows(currentRows(i)));
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

    RF_Pred_valid(currentRows, :) = blockPred;
    RF_Val_Pred_valid(currentRows, :) = blockValPred;
    RF_NRMSE_valid(currentRows, :) = blockNRMSE;
    RF_TDA_valid(currentRows, :) = blockTDA;
    RF_best_p_valid(currentRows, :) = blockBestP;
    RF_Val_NRMSE_best_valid(currentRows, :) = blockValBest;
    RF_RMSE_valid(currentRows, :) = blockRMSE;
    RF_TestRange_valid(currentRows, :) = blockTestRange;
    RF_TestSuccess_valid(currentRows, :) = blockSuccess;
    RF_FailureCode_valid(currentRows, :) = blockFailureCode;
    RF_FailureID_valid(currentRows, :) = blockFailureID;
    blockDone(block) = true;

    checkpointInfo = struct('version', checkpointVersion, ...
        'runSignature', runSignature);
    checkpointTemp = fullfile(outFolder, [checkpointFileName, '.tmp.mat']);
    save(checkpointTemp, ...
        'RF_Pred_valid', 'RF_Val_Pred_valid', ...
        'RF_NRMSE_valid', 'RF_TDA_valid', 'RF_best_p_valid', ...
        'RF_Val_NRMSE_best_valid', 'RF_RMSE_valid', 'RF_TestRange_valid', ...
        'RF_TestSuccess_valid', 'RF_FailureCode_valid', 'RF_FailureID_valid', ...
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
fprintf('\nRF computation finished.\n');
fprintf('Computation time this run: %.2f minutes (%.2f hours)\n', ...
    elapsedComputeThisRun / 60, elapsedComputeThisRun / 3600);
fprintf('Cumulative block computation time: %.2f minutes (%.2f hours)\n', ...
    cumulativeComputeTimeSeconds / 60, cumulativeComputeTimeSeconds / 3600);
if useParallel
    delete(Y_valid_const);
end

%% 3.5 Map valid-grid results to 1-degree fields
RF_Pred = mapValidToGrid3D(RF_Pred_valid, validLinearIndex, nLat, nLon);
RF_Val_Pred = mapValidToGrid3D(RF_Val_Pred_valid, validLinearIndex, nLat, nLon);
RF_NRMSE = mapValidToGrid3D(RF_NRMSE_valid, validLinearIndex, nLat, nLon);
RF_TDA = mapValidToGrid3D(RF_TDA_valid, validLinearIndex, nLat, nLon);
RF_best_p = mapValidToGrid3D(RF_best_p_valid, validLinearIndex, nLat, nLon);
RF_Val_NRMSE_best = mapValidToGrid3D(RF_Val_NRMSE_best_valid, ...
    validLinearIndex, nLat, nLon);
RF_RMSE = mapValidToGrid3D(RF_RMSE_valid, validLinearIndex, nLat, nLon);
RF_TestRange = mapValidToGrid3D(RF_TestRange_valid, validLinearIndex, nLat, nLon);
RF_TestSuccess = mapValidToGrid3D(double(RF_TestSuccess_valid), ...
    validLinearIndex, nLat, nLon);
RF_FailureCode = mapValidToGrid3D(double(RF_FailureCode_valid), ...
    validLinearIndex, nLat, nLon);

RF_NRMSE_mean = mean(RF_NRMSE, 3, 'omitnan');
RF_TDA_mean = mean(RF_TDA, 3, 'omitnan');
RF_Val_NRMSE_best_mean = mean(RF_Val_NRMSE_best, 3, 'omitnan');

CompleteThreeYears_valid = all(RF_TestSuccess_valid, 2);
RF_CompleteThreeYears = false(nLat, nLon);
RF_CompleteThreeYears(validLinearIndex) = CompleteThreeYears_valid;
RF_NRMSE_mean_complete = RF_NRMSE_mean;
RF_NRMSE_mean_complete(~RF_CompleteThreeYears) = NaN;

fprintf('\nRF global results:\n');
fprintf('Successful grid-year forecasts: %d / %d.\n', ...
    nnz(RF_TestSuccess_valid), numel(RF_TestSuccess_valid));
fprintf('Complete three-year grid cells: %d / %d.\n', ...
    nnz(CompleteThreeYears_valid), nValid);
for r = 1:nRound
    fprintf('  %d success: %d / %d\n', testYears(r), ...
        nnz(RF_TestSuccess_valid(:, r)), nValid);
    for code = 1:4
        nFailure = nnz(RF_FailureCode_valid(:, r) == code);
        if nFailure > 0
            fprintf('    Failure code %d: %d\n', code, nFailure);
        end
    end
end
fprintf('RF_NRMSE_mean (finite values): %.4f\n', ...
    mean(RF_NRMSE_mean(CommonValid), 'omitnan'));
fprintf('RF_TDA_mean (finite values): %.2f%%\n', ...
    mean(RF_TDA_mean(CommonValid), 'omitnan'));

%% 3.6 Save final output
Meta = struct();
Meta.model = 'RF';
Meta.experimentName = experimentName;
Meta.input = 'Monthly dGWS from the common GLDAS valid-grid database';
Meta.period = '2004-01 to 2024-12';
Meta.predictionPeriod = '2022-01 to 2024-12';
Meta.validationPeriod = '2021-01 to 2023-12';
Meta.arrayConvention = 'Maps: [latitude, longitude, time]; valid series: [grid, time]';
Meta.rollingStrategy = ['Train 2004 through test year minus two; ' ...
    'validate in the following year; test in the subsequent year'];
Meta.p_list = p_list;
Meta.rfParam = rfParam;
Meta.numRepeat = numRepeat;
Meta.randomSeed = randomSeed;
Meta.standardization = 'None; regression trees use original dGWS lagged values';
Meta.forecasting = 'Twelve direct RF ensembles, one for each lead month';
Meta.finalPrediction = ['The fitted forests selected on the validation year ' ...
    'are reused for test forecasting without refitting. Observed ' ...
    'validation-year values update only the input history window'];
Meta.repeatStrategy = ['Validation selects p by the mean of finite replicate ' ...
    'NRMSE values. Validation and test predictions average the fitted ' ...
    'forests across successful replicates. For numRepeat > 1, NRMSE ' ...
    'of the averaged validation prediction is not the selection score'];
Meta.randomization = ['Each forest fit receives a deterministic seed based ' ...
    'on experiment, grid row, test round, history length, replicate ' ...
    'and horizon. Tree predictor selection is reproducible. ' ...
    'No additional random fit is made for testing'];
Meta.NRMSE = 'Annual RMSE divided by observed test-year dGWS range';
Meta.TDA = ['Percentage of matching signs of successive month-to-month ' ...
    'differences in dGWS within a test year (11 comparisons)'];
Meta.failureCodes = ['0 = success; 1 = no finite validation candidate; ' ...
    '2 = test forecast error; 3 = invalid test predictions; ' ...
    '4 = nonfinite test NRMSE'];
Meta.RF_FailureID = 'MATLAB exception identifiers for failed model fits or forecasts';
Meta.meanMetricRule = 'Three-year mean computed with omitnan';
Meta.completeMetricRule = ['RF_NRMSE_mean_complete requires successful forecasts ' ...
    'in all three test years'];
Meta.nValid = nValid;
Meta.block_num = block_num;
Meta.numWorkers = numWorkers;
Meta.usedParallel = useParallel;
Meta.runSignature = runSignature;
Meta.computeTimeMinutes = cumulativeComputeTimeSeconds / 60;
Meta.computeTimeThisRunMinutes = elapsedComputeThisRun / 60;
Meta.scriptName = 'RF_global_forecast.m';
Meta.pilotTag = cfg.tag;
Meta.pilotOriginalValidRows = originalValidRows;
Meta.pilotN = nValid;
Meta.pilotMode = 'complete-revalidation';

finalTemp = fullfile(outFolder, [outFileName, '.tmp.mat']);
fprintf('\nSaving final MAT-file: %s\n', outFile);
save(finalTemp, ...
    'RF_Pred', 'RF_Val_Pred', ...
    'RF_NRMSE', 'RF_NRMSE_mean', 'RF_NRMSE_mean_complete', ...
    'RF_TDA', 'RF_TDA_mean', 'RF_best_p', ...
    'RF_Val_NRMSE_best', 'RF_Val_NRMSE_best_mean', ...
    'RF_RMSE', 'RF_TestRange', 'RF_TestSuccess', ...
    'RF_FailureCode', 'RF_CompleteThreeYears', ...
    'RF_Pred_valid', 'RF_Val_Pred_valid', ...
    'RF_NRMSE_valid', 'RF_TDA_valid', 'RF_best_p_valid', ...
    'RF_Val_NRMSE_best_valid', 'RF_RMSE_valid', ...
    'RF_TestRange_valid', 'RF_TestSuccess_valid', ...
    'RF_FailureCode_valid', 'RF_FailureID_valid', ...
    'p_list', 'numRepeat', 'testYears', 'predYM', 'valYM', ...
    'rfParam', 'randomSeed', 'experimentName', ...
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
    figure('Name', 'RF mean NRMSE');
    p = imagesc(lon, lat, RF_NRMSE_mean);
    set(p, 'AlphaData', double(isfinite(RF_NRMSE_mean)));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(jet);
    clim(plotNRMSELimits);
    cb = colorbar;
    ylabel(cb, 'NRMSE');
    title('RF mean NRMSE, test years 2022-2024');
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

elapsedTotalSeconds = toc(totalTic);
fprintf('\nRF processing completed.\n');
fprintf('Total time since computation started (including output): %.2f minutes (%.2f hours)\n', ...
    elapsedTotalSeconds / 60, elapsedTotalSeconds / 3600);

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

end  % End of SKLGB_pilot_RF main function

function [pred36, valPred36, nrmse3, tda3, bestP3, valBest3, ...
    rmse3, testRange3, success3, failureCode3, failureID3] = ...
    runRFOnePixelAllRounds(Y, p_list, numRepeat, testYears, horizon, ...
    rfParam, randomSeed, pointID)

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

    if numel(Y) ~= 252 || any(~isfinite(Y))
        error('RF input must be a finite, 252-month time series.');
    end

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
        bestValScore = Inf;
        selectedPacks = cell(numRepeat, 1);
        selectedValPred = NaN(horizon, 1);
        bestP = NaN;
        firstValidationError = '';

        for k = 1:numel(p_list)
            p = p_list(k);
            if numel(Y_train) <= p + horizon - 1
                continue;
            end

            valScores = NaN(numRepeat, 1);
            candidatePacks = cell(numRepeat, 1);
            candidateValPred = NaN(horizon, numRepeat);
            initWindowVal = Y_train(end-p+1:end);

            for rep = 1:numRepeat
                try
                    seedVal = makeSeed(randomSeed, pointID, r, k, rep);
                    pack = fitRFHorizonSpecificFromSeries( ...
                        Y_train, p, horizon, rfParam, seedVal);
                    predVal = predictRFHorizonSpecific(pack, initWindowVal);
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

            % Each replicate has the same training set and candidate p.
            % The selection criterion is mean individual validation NRMSE.
            validRep = isfinite(valScores);
            if ~any(validRep)
                continue;
            end
            candidateScore = mean(valScores(validRep));
            if isfinite(candidateScore) && candidateScore < bestValScore
                bestValScore = candidateScore;
                bestP = p;
                selectedPacks = candidatePacks;
                selectedValPred = mean(candidateValPred(:, validRep), 2);
            end
        end

        if ~isfinite(bestValScore)
            failureCode3(r) = 1;
            failureID3{r} = firstValidationError;
            continue;
        end

        bestP3(r) = bestP;
        valBest3(r) = bestValScore;
        valPred36(predIdx) = selectedValPred(:);

        % Test from the next December using the same fitted forests.
        initWindowTest = Y(trainEndIdx + horizon - bestP + 1: ...
            trainEndIdx + horizon);
        testPredRepeats = NaN(horizon, numRepeat);
        firstTestError = '';
        for rep = 1:numRepeat
            pack = selectedPacks{rep};
            if isempty(pack)
                continue;
            end
            try
                testPredRepeats(:, rep) = predictRFHorizonSpecific( ...
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

function rfPack = fitRFHorizonSpecificFromSeries( ...
    Y_series, p, horizon, rfParam, baseSeed)

    Y_series = double(Y_series(:));
    modelCell = cell(horizon, 1);
    numVariablesToSample = getNumVariablesToSample(p, rfParam);
    treeTemplate = templateTree( ...
        'MinLeafSize', rfParam.MinLeafSize, ...
        'NumVariablesToSample', numVariablesToSample, ...
        'Reproducible', true);

    for h = 1:horizon
        [X, Y_target] = buildLaggedSamplesDirect(Y_series, p, h);
        seedH = mod(baseSeed + h, 2^32 - 1);
        if seedH == 0
            seedH = 1;
        end
        rng(seedH, 'twister');
        modelCell{h} = fitrensemble(X, Y_target, ...
            'Method', rfParam.Method, ...
            'NumLearningCycles', rfParam.NumLearningCycles, ...
            'Learners', treeTemplate);
    end

    rfPack = struct();
    rfPack.modelCell = modelCell;
    rfPack.p = p;
    rfPack.horizon = horizon;
end

function [X, Y_target] = buildLaggedSamplesDirect(Y_series, p, h)
    Y_series = double(Y_series(:));
    nSample = numel(Y_series) - p - h + 1;
    if nSample <= 0
        error('No direct samples can be constructed.');
    end
    X = zeros(nSample, p);
    Y_target = zeros(nSample, 1);
    for j = 1:nSample
        X(j, :) = Y_series(j:j+p-1)';
        Y_target(j) = Y_series(j+p+h-1);
    end
end

function Y_pred = predictRFHorizonSpecific(rfPack, initialWindow)
    p = rfPack.p;
    horizon = rfPack.horizon;
    initialWindow = double(initialWindow(:)');
    if numel(initialWindow) ~= p
        error('Initial history window length does not match p.');
    end
    Y_pred = NaN(horizon, 1);
    for h = 1:horizon
        prediction = predict(rfPack.modelCell{h}, initialWindow);
        Y_pred(h) = double(prediction(1));
    end
end

function nVars = getNumVariablesToSample(p, rfParam)
    switch rfParam.NumVariablesToSampleMode
        case 'p_over_3'
            nVars = max(1, floor(p / 3));
        case 'sqrt'
            nVars = max(1, floor(sqrt(p)));
        case 'all'
            nVars = p;
        otherwise
            error('Unknown NumVariablesToSampleMode.');
    end
    nVars = min(nVars, p);
end

function seedVal = makeSeed(randomSeed, pointID, roundID, pID, repeatID)
    seedRaw = double(randomSeed) ...
        + double(pointID) * 100 ...
        + double(roundID) * 10000 ...
        + double(pID) * 1000 ...
        + double(repeatID) * 10;
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
    p_list, testYears, horizon, numRepeat, rfParam, randomSeed, blockNum, version)
    md = java.security.MessageDigest.getInstance('SHA-256');
    updateDigest(md, 'RF_global_forecast');
    updateDigest(md, double(version));
    updateDigest(md, data);
    updateDigest(md, index);
    updateDigest(md, lat);
    updateDigest(md, lon);
    updateDigest(md, ym);
    updateDigest(md, double(p_list));
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
        'RF_Pred_valid', 'RF_Val_Pred_valid', ...
        'RF_NRMSE_valid', 'RF_TDA_valid', ...
        'RF_best_p_valid', 'RF_Val_NRMSE_best_valid', ...
        'RF_RMSE_valid', 'RF_TestRange_valid', 'RF_TestSuccess_valid', ...
        'RF_FailureCode_valid', 'RF_FailureID_valid', ...
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
    matrix3 = {'RF_NRMSE_valid', 'RF_TDA_valid', ...
        'RF_best_p_valid', 'RF_Val_NRMSE_best_valid', ...
        'RF_RMSE_valid', 'RF_TestRange_valid', 'RF_TestSuccess_valid', ...
        'RF_FailureCode_valid', 'RF_FailureID_valid'};
    for k = 1:numel(matrix3)
        if ~isequal(size(C.(matrix3{k})), [nValid, nRound])
            error('Checkpoint contains an invalid %s array.', matrix3{k});
        end
    end
    if ~isequal(size(C.RF_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.RF_Val_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.blockDone), [totalBlocks, 1])
        error('Checkpoint dimensions do not match the current run.');
    end
end
