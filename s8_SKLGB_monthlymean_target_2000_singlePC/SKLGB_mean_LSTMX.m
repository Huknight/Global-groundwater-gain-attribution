function SKLGB_mean_LSTMX(cfg)
% PILOT WRAPPER: adapted from LSTMX_global_forecast.m.
% Only file paths, pilot seed IDs, run settings and the two tested
% hyperparameters are changed; the original model code follows.
% Each configuration searches the FULL original validation p/x grid.

%% 1. CODE DESCRIPTION
% LSTMX - Global multivariate groundwater forecasting
% Purpose
%   Forecast monthly groundwater-storage change (dGWS) at each valid
%   1-degree cell using a sequence-to-vector LSTM with historical inputs.
%
% Input
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   dGWS_valid, Prec_valid, Evap_valid, Runoff_valid, dSM2m_valid,
%   SurfT_valid: [grid, month], January 2004 through December 2024.
%
% Rolling evaluation
%   Test 2022: train 2004-2020, validate 2021.
%   Test 2023: train 2004-2021, validate 2022.
%   Test 2024: train 2004-2022, validate 2023.
%   Forecasts use observations available through the preceding December;
%   no future exogenous observations enter validation or test forecasts.
%
% Model
%   Direct 12-month sequence-to-vector regression: sequenceInputLayer(6),
%   LSTM with 32 hidden units and OutputMode='last', fullyConnectedLayer(12),
%   regressionLayer. Common historical window p in [1, 6, 12, 24, 36, 48, 60]
%   for dGWS and all five exogenous variables.
%   Training-only normalization is applied by predictor variable, across
%   training windows and time positions, and by output horizon for dGWS.
%   Validation-year observations control early stopping and select p by NRMSE.
%   Training returns the fitted network with the lowest validation loss.
%   The selected network is reused for testing without refitting; only
%   the forecast-origin historical input is updated.
%   When numRepeat > 1, independent selected networks are averaged.
%
% Evaluation
%   NRMSE = RMSE / (max(test dGWS) - min(test dGWS)).
%   TDA compares signs of successive changes in monthly dGWS,
%   using eleven comparisons within each test year.
%
% Output
%   LSTMX_Pred, LSTMX_Val_Pred: [lat, lon, 36], chronological months.
%   LSTMX_NRMSE, LSTMX_TDA, LSTMX_best_p: [lat, lon, 3].
%   LSTMX_Val_NRMSE_best, LSTMX_RMSE, LSTMX_TestRange: [lat, lon, 3].
%   LSTMX_TestSuccess, LSTMX_FailureCode: [lat, lon, 3].
%   Additional valid-grid outputs, coordinates and metadata are saved.
%
%% 2. CENTRALIZED PARAMETERS
% Pilot wrapper: do not clear the input configuration.
projectRoot = cfg.projectRoot;
dataFile = cfg.pilotMasterFile;
outFolder = fullfile(cfg.outputRoot, 'LSTMX');

experimentName = 'original';  % 'original', 'uc1', or 'uc2'
switch lower(experimentName)
    case 'original'
        randomSeed = 20260605;
        outFileName = 'LSTMX_monthmean_2022_2024.mat';
        checkpointFileName = 'LSTMX_monthmean_2022_2024_checkpoint.mat';
    case 'uc1'
        randomSeed = 20260719;
        outFileName = 'LSTMX_multivariate_global_200401_202412_uc1.mat';
        checkpointFileName = 'LSTMX_multivariate_global_uc1_checkpoint.mat';
    case 'uc2'
        randomSeed = 20260721;
        outFileName = 'LSTMX_multivariate_global_200401_202412_uc2.mat';
        checkpointFileName = 'LSTMX_multivariate_global_uc2_checkpoint.mat';
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
numRepeat = 1;
inputVarNames = {'dGWS', 'Prec', 'Evap', 'Runoff', 'dSM2m', 'SurfT'};

lstmParam = struct();
lstmParam.numPredict = 12;
lstmParam.numFeatures = numel(inputVarNames);
lstmParam.numHiddenUnits = 32;
lstmParam.initialLearnRate = 0.01;
lstmParam.miniBatchSize = 8;
lstmParam.maxEpochs = 500;
lstmParam.validationPatience = 10;
lstmParam.outputNetwork = 'best-validation';
lstmParam.gradientThreshold = 1;
lstmParam.optimizer = 'adam';
lstmParam.executionEnvironment = 'cpu';
lstmParam.shuffle = 'never';
forecastHorizon = lstmParam.numPredict;

makeDiagnosticPlot = false;
plotNRMSELimits = [0, 0.4];
checkpointVersion = 3;

%% 3. MAIN EXECUTION

% 3.1 Check settings and read the master database
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
if lstmParam.numPredict ~= 12 || lstmParam.numFeatures ~= numel(inputVarNames) || ...
        ~isscalar(lstmParam.numHiddenUnits) || lstmParam.numHiddenUnits < 1 || ...
        lstmParam.numHiddenUnits ~= floor(lstmParam.numHiddenUnits) || ...
        ~isscalar(lstmParam.initialLearnRate) || ...
        lstmParam.initialLearnRate <= 0 || ...
        ~isscalar(lstmParam.miniBatchSize) || lstmParam.miniBatchSize < 1 || ...
        lstmParam.miniBatchSize ~= floor(lstmParam.miniBatchSize) || ...
        ~isscalar(lstmParam.maxEpochs) || lstmParam.maxEpochs < 1 || ...
        lstmParam.maxEpochs ~= floor(lstmParam.maxEpochs) || ...
        ~isscalar(lstmParam.validationPatience) || ...
        lstmParam.validationPatience < 1 || ...
        lstmParam.validationPatience ~= floor(lstmParam.validationPatience) || ...
        ~strcmpi(lstmParam.outputNetwork, 'best-validation') || ...
        ~isscalar(lstmParam.gradientThreshold) || ...
        lstmParam.gradientThreshold <= 0 || ...
        ~strcmpi(lstmParam.optimizer, 'adam') || ...
        ~strcmpi(lstmParam.executionEnvironment, 'cpu') || ...
        ~strcmpi(lstmParam.shuffle, 'never')
    error('Invalid LSTMX training settings.');
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

requiredVars = {'dGWS_valid', ...
    'Prec_valid', 'Evap_valid', 'Runoff_valid', 'dSM2m_valid', 'SurfT_valid', ...
    'CommonMask', 'CommonValid', ...
    'validLinearIndex', 'validLat', 'validLon', 'lat', 'lon', 'ym'};
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
dataVars = {dGWS_valid, Prec_valid, Evap_valid, Runoff_valid, ...
    dSM2m_valid, SurfT_valid};
for v = 1:numel(dataVars)
    V = dataVars{v};
    if ~isequal(size(V), [nValid, nTime]) || nTime ~= 252 || ...
            any(~isfinite(V(:)))
        error('%s must contain finite [grid, 252] time series.', ...
            inputVarNames{v});
    end
end
clear V dataVars
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
fprintf('LSTMX experiment: %s, random seed %d.\n', experimentName, randomSeed);
fprintf('Inputs: %s.\n', strjoin(inputVarNames, ', '));
fprintf('Grid cells: %d; monthly observations: %d; test rounds: %d.\n', ...
    nValid, nTime, nRound);
fprintf('History lengths: %s\n', mat2str(p_list));
fprintf('LSTMX hidden units: %d; epochs: %d; patience: %d; batch size: %d.\n', ...
    lstmParam.numHiddenUnits, lstmParam.maxEpochs, ...
    lstmParam.validationPatience, lstmParam.miniBatchSize);
fprintf('Optimizer: %s; initial learning rate: %.4g; execution: %s.\n', ...
    lstmParam.optimizer, lstmParam.initialLearnRate, ...
    lstmParam.executionEnvironment);
fprintf('Output network: %s (minimum validation loss).\n', ...
    lstmParam.outputNetwork);
fprintf('Replicates per candidate: %d.\n', numRepeat);
fprintf('Blocks: %d; approximately %d grid cells per block.\n', ...
    total_blocks, block_size);

% 3.2 Prepare checkpoint and result arrays
runSignature = computeRunSignature(dGWS_valid, Prec_valid, Evap_valid, ...
    Runoff_valid, dSM2m_valid, SurfT_valid, validLinearIndex, ...
    lat, lon, ym, p_list, testYears, forecastHorizon, ...
    numRepeat, lstmParam, randomSeed, block_num, checkpointVersion);

if isfile(checkpointFile) && ~resumeFromCheckpoint
    delete(checkpointFile);
end

if resumeFromCheckpoint && isfile(checkpointFile)
    fprintf('Reading checkpoint: %s\n', checkpointFile);
    C = load(checkpointFile);
    validateCheckpoint(C, runSignature, nValid, nRound, nPredMonth, total_blocks);

    LSTMX_Pred_valid = C.LSTMX_Pred_valid;
    LSTMX_Val_Pred_valid = C.LSTMX_Val_Pred_valid;
    LSTMX_NRMSE_valid = C.LSTMX_NRMSE_valid;
    LSTMX_TDA_valid = C.LSTMX_TDA_valid;
    LSTMX_best_p_valid = C.LSTMX_best_p_valid;
    LSTMX_Val_NRMSE_best_valid = C.LSTMX_Val_NRMSE_best_valid;
    LSTMX_RMSE_valid = C.LSTMX_RMSE_valid;
    LSTMX_TestRange_valid = C.LSTMX_TestRange_valid;
    LSTMX_TestSuccess_valid = C.LSTMX_TestSuccess_valid;
    LSTMX_FailureCode_valid = C.LSTMX_FailureCode_valid;
    LSTMX_FailureID_valid = C.LSTMX_FailureID_valid;
    blockDone = C.blockDone;
    cumulativeComputeTimeSeconds = C.cumulativeComputeTimeSeconds;
    clear C
else
    LSTMX_Pred_valid = NaN(nValid, nPredMonth);
    LSTMX_Val_Pred_valid = NaN(nValid, nPredMonth);
    LSTMX_NRMSE_valid = NaN(nValid, nRound);
    LSTMX_TDA_valid = NaN(nValid, nRound);
    LSTMX_best_p_valid = NaN(nValid, nRound);
    LSTMX_Val_NRMSE_best_valid = NaN(nValid, nRound);
    LSTMX_RMSE_valid = NaN(nValid, nRound);
    LSTMX_TestRange_valid = NaN(nValid, nRound);
    LSTMX_TestSuccess_valid = false(nValid, nRound);
    LSTMX_FailureCode_valid = uint8(zeros(nValid, nRound));
    LSTMX_FailureID_valid = cell(nValid, nRound);
    blockDone = false(total_blocks, 1);
    cumulativeComputeTimeSeconds = 0;
end

% 3.3 Start parallel pool
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
    dGWS_const = parallel.pool.Constant(dGWS_valid);
    Prec_const = parallel.pool.Constant(Prec_valid);
    Evap_const = parallel.pool.Constant(Evap_valid);
    Runoff_const = parallel.pool.Constant(Runoff_valid);
    dSM2m_const = parallel.pool.Constant(dSM2m_valid);
    SurfT_const = parallel.pool.Constant(SurfT_valid);
end

% 3.4 Run global LSTMX forecasts
fprintf('\nStarting LSTMX computation...\n');
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
            X_pixel = [
                dGWS_const.Value(currentRows(i), :)', ...
                Prec_const.Value(currentRows(i), :)', ...
                Evap_const.Value(currentRows(i), :)', ...
                Runoff_const.Value(currentRows(i), :)', ...
                dSM2m_const.Value(currentRows(i), :)', ...
                SurfT_const.Value(currentRows(i), :)'];
            Y_pixel = X_pixel(:, 1);
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runLSTMXOnePixelAllRounds(X_pixel, Y_pixel, ...
                p_list, numRepeat, testYears, forecastHorizon, lstmParam, ...
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
            X_pixel = [
                dGWS_valid(currentRows(i), :)', ...
                Prec_valid(currentRows(i), :)', ...
                Evap_valid(currentRows(i), :)', ...
                Runoff_valid(currentRows(i), :)', ...
                dSM2m_valid(currentRows(i), :)', ...
                SurfT_valid(currentRows(i), :)'];
            Y_pixel = X_pixel(:, 1);
            [p36, v36, n3, t3, best3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runLSTMXOnePixelAllRounds(X_pixel, Y_pixel, ...
                p_list, numRepeat, testYears, forecastHorizon, lstmParam, ...
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

    LSTMX_Pred_valid(currentRows, :) = blockPred;
    LSTMX_Val_Pred_valid(currentRows, :) = blockValPred;
    LSTMX_NRMSE_valid(currentRows, :) = blockNRMSE;
    LSTMX_TDA_valid(currentRows, :) = blockTDA;
    LSTMX_best_p_valid(currentRows, :) = blockBestP;
    LSTMX_Val_NRMSE_best_valid(currentRows, :) = blockValBest;
    LSTMX_RMSE_valid(currentRows, :) = blockRMSE;
    LSTMX_TestRange_valid(currentRows, :) = blockTestRange;
    LSTMX_TestSuccess_valid(currentRows, :) = blockSuccess;
    LSTMX_FailureCode_valid(currentRows, :) = blockFailureCode;
    LSTMX_FailureID_valid(currentRows, :) = blockFailureID;
    blockDone(block) = true;

    checkpointInfo = struct('version', checkpointVersion, ...
        'runSignature', runSignature);
    checkpointTemp = fullfile(outFolder, [checkpointFileName, '.tmp.mat']);
    save(checkpointTemp, ...
        'LSTMX_Pred_valid', 'LSTMX_Val_Pred_valid', ...
        'LSTMX_NRMSE_valid', 'LSTMX_TDA_valid', 'LSTMX_best_p_valid', ...
        'LSTMX_Val_NRMSE_best_valid', 'LSTMX_RMSE_valid', 'LSTMX_TestRange_valid', ...
        'LSTMX_TestSuccess_valid', 'LSTMX_FailureCode_valid', 'LSTMX_FailureID_valid', ...
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
fprintf('\nLSTMX computation finished.\n');
fprintf('Computation time this run: %.2f minutes (%.2f hours)\n', ...
    elapsedComputeThisRun / 60, elapsedComputeThisRun / 3600);
fprintf('Cumulative block computation time: %.2f minutes (%.2f hours)\n', ...
    cumulativeComputeTimeSeconds / 60, cumulativeComputeTimeSeconds / 3600);
if useParallel
    delete(dGWS_const);
    delete(Prec_const);
    delete(Evap_const);
    delete(Runoff_const);
    delete(dSM2m_const);
    delete(SurfT_const);
end

% 3.5 Map valid-grid results to 1-degree fields
LSTMX_Pred = mapValidToGrid3D(LSTMX_Pred_valid, validLinearIndex, nLat, nLon);
LSTMX_Val_Pred = mapValidToGrid3D(LSTMX_Val_Pred_valid, validLinearIndex, nLat, nLon);
LSTMX_NRMSE = mapValidToGrid3D(LSTMX_NRMSE_valid, validLinearIndex, nLat, nLon);
LSTMX_TDA = mapValidToGrid3D(LSTMX_TDA_valid, validLinearIndex, nLat, nLon);
LSTMX_best_p = mapValidToGrid3D(LSTMX_best_p_valid, validLinearIndex, nLat, nLon);
LSTMX_Val_NRMSE_best = mapValidToGrid3D(LSTMX_Val_NRMSE_best_valid, ...
    validLinearIndex, nLat, nLon);
LSTMX_RMSE = mapValidToGrid3D(LSTMX_RMSE_valid, validLinearIndex, nLat, nLon);
LSTMX_TestRange = mapValidToGrid3D(LSTMX_TestRange_valid, validLinearIndex, nLat, nLon);
LSTMX_TestSuccess = mapValidToGrid3D(double(LSTMX_TestSuccess_valid), ...
    validLinearIndex, nLat, nLon);
LSTMX_FailureCode = mapValidToGrid3D(double(LSTMX_FailureCode_valid), ...
    validLinearIndex, nLat, nLon);

LSTMX_NRMSE_mean = mean(LSTMX_NRMSE, 3, 'omitnan');
LSTMX_TDA_mean = mean(LSTMX_TDA, 3, 'omitnan');
LSTMX_Val_NRMSE_best_mean = mean(LSTMX_Val_NRMSE_best, 3, 'omitnan');

CompleteThreeYears_valid = all(LSTMX_TestSuccess_valid, 2);
LSTMX_CompleteThreeYears = false(nLat, nLon);
LSTMX_CompleteThreeYears(validLinearIndex) = CompleteThreeYears_valid;
LSTMX_NRMSE_mean_complete = LSTMX_NRMSE_mean;
LSTMX_NRMSE_mean_complete(~LSTMX_CompleteThreeYears) = NaN;

fprintf('\nLSTMX global results:\n');
fprintf('Successful grid-year forecasts: %d / %d.\n', ...
    nnz(LSTMX_TestSuccess_valid), numel(LSTMX_TestSuccess_valid));
fprintf('Complete three-year grid cells: %d / %d.\n', ...
    nnz(CompleteThreeYears_valid), nValid);
for r = 1:nRound
    fprintf('  %d success: %d / %d\n', testYears(r), ...
        nnz(LSTMX_TestSuccess_valid(:, r)), nValid);
    for code = 1:4
        nFailure = nnz(LSTMX_FailureCode_valid(:, r) == code);
        if nFailure > 0
            fprintf('    Failure code %d: %d\n', code, nFailure);
        end
    end
end
fprintf('LSTMX_NRMSE_mean (finite values): %.4f\n', ...
    mean(LSTMX_NRMSE_mean(CommonValid), 'omitnan'));
fprintf('LSTMX_TDA_mean (finite values): %.2f%%\n', ...
    mean(LSTMX_TDA_mean(CommonValid), 'omitnan'));

% 3.6 Save final output
Meta = struct();
Meta.targetDefinition='difference of consecutive 1-degree monthly-mean GWS states';
Meta.model = 'LSTMX';
Meta.experimentName = experimentName;
Meta.input = ['Historical dGWS, precipitation, evapotranspiration, runoff, ' ...
    'soil-moisture change and surface temperature'];
Meta.target = 'Future monthly dGWS';
Meta.period = '2004-01 to 2024-12';
Meta.predictionPeriod = '2022-01 to 2024-12';
Meta.validationPeriod = '2021-01 to 2023-12';
Meta.arrayConvention = 'Maps: [latitude, longitude, time]; valid series: [grid, time]';
Meta.rollingStrategy = ['Train 2004 through test year minus two; ' ...
    'validate in the following year; test in the subsequent year'];
Meta.p_list = p_list;
Meta.lstmParam = lstmParam;
Meta.numRepeat = numRepeat;
Meta.randomSeed = randomSeed;
Meta.inputVarNames = inputVarNames;
Meta.standardization = ['Training-only feature-wise normalization across training ' ...
    'sequences and time positions, and horizon-wise response normalization'];
Meta.architecture = ['sequenceInputLayer(6), lstmLayer(32, OutputMode=last), ' ...
    'fullyConnectedLayer(12), regressionLayer'];
Meta.forecasting = 'One direct six-input LSTM predicts the next 12 months';
Meta.noFutureExogenous = 'Only predictors observed at or before the forecast origin are used';
Meta.earlyStopping = ['Validation-year targets are used for early stopping; ' ...
    'they are not included in gradient-based fitting samples'];
Meta.outputNetwork = ['The network with the lowest validation loss is returned ' ...
    'and reused for testing without refitting'];
Meta.finalPrediction = ['The fitted LSTMX selected on the validation year ' ...
    'is reused for test forecasting without refitting. Observed ' ...
    'validation-year values update only the input history window'];
Meta.repeatStrategy = ['Validation selects p by mean finite replicate NRMSE. ' ...
    'For numRepeat greater than one, validation and test predictions ' ...
    'average the corresponding successful fitted networks'];
Meta.randomization = ['Network training uses deterministic seeds based on ' ...
    'experiment, grid row, test round, history length, and replicate. ' ...
    'No additional random initialization is made for testing'];
Meta.NRMSE = 'Annual RMSE divided by observed test-year dGWS range';
Meta.TDA = ['Percentage of matching signs of successive month-to-month ' ...
    'differences in dGWS within a test year (11 comparisons)'];
Meta.failureCodes = ['0 = success; 1 = no finite validation candidate; ' ...
    '2 = test forecast error; 3 = invalid test predictions; ' ...
    '4 = nonfinite test NRMSE'];
Meta.LSTMX_FailureID = 'MATLAB exception identifiers for failed model fits or forecasts';
Meta.meanMetricRule = 'Three-year mean computed with omitnan';
Meta.completeMetricRule = ['LSTMX_NRMSE_mean_complete requires successful forecasts ' ...
    'in all three test years'];
Meta.nValid = nValid;
Meta.block_num = block_num;
Meta.numWorkers = numWorkers;
Meta.usedParallel = useParallel;
Meta.runSignature = runSignature;
Meta.computeTimeMinutes = cumulativeComputeTimeSeconds / 60;
Meta.computeTimeThisRunMinutes = elapsedComputeThisRun / 60;
Meta.scriptName = 'SKLGB_mean_LSTMX.m';
Meta.pilotTag = 'target_monthmean_2022_2024';
Meta.pilotOriginalValidRows = originalValidRows;
Meta.pilotN = nValid;
Meta.pilotMode = 'monthlymean-target-2022-2024-original-hyperparameters';

finalTemp = fullfile(outFolder, [outFileName, '.tmp.mat']);
fprintf('\nSaving final MAT-file: %s\n', outFile);
save(finalTemp, ...
    'LSTMX_Pred', 'LSTMX_Val_Pred', ...
    'LSTMX_NRMSE', 'LSTMX_NRMSE_mean', 'LSTMX_NRMSE_mean_complete', ...
    'LSTMX_TDA', 'LSTMX_TDA_mean', 'LSTMX_best_p', ...
    'LSTMX_Val_NRMSE_best', 'LSTMX_Val_NRMSE_best_mean', ...
    'LSTMX_RMSE', 'LSTMX_TestRange', 'LSTMX_TestSuccess', ...
    'LSTMX_FailureCode', 'LSTMX_CompleteThreeYears', ...
    'LSTMX_Pred_valid', 'LSTMX_Val_Pred_valid', ...
    'LSTMX_NRMSE_valid', 'LSTMX_TDA_valid', 'LSTMX_best_p_valid', ...
    'LSTMX_Val_NRMSE_best_valid', 'LSTMX_RMSE_valid', ...
    'LSTMX_TestRange_valid', 'LSTMX_TestSuccess_valid', ...
    'LSTMX_FailureCode_valid', 'LSTMX_FailureID_valid', ...
    'p_list', 'numRepeat', 'testYears', 'predYM', 'valYM', ...
    'lstmParam', 'randomSeed', 'experimentName', 'inputVarNames', ...
    'CommonMask', 'CommonValid', 'validLinearIndex', 'nValid', ...
    'validLat', 'validLon', 'lat', 'lon', 'ym', 'yearVec', 'monthVec', ...
    'Meta', '-v7.3');
[moved, msg] = movefile(finalTemp, outFile, 'f');
if ~moved
    error('Final MAT-file replacement failed: %s', msg);
end
fprintf('Saved successfully.\n');

% 3.7 Diagnostic map
if makeDiagnosticPlot
    figure('Name', 'LSTMX mean NRMSE');
    p = imagesc(lon, lat, LSTMX_NRMSE_mean);
    set(p, 'AlphaData', double(isfinite(LSTMX_NRMSE_mean)));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(jet);
    clim(plotNRMSELimits);
    cb = colorbar;
    ylabel(cb, 'NRMSE');
    title('LSTMX mean NRMSE, test years 2022-2024');
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

elapsedTotalSeconds = toc(totalTic);
fprintf('\nLSTMX processing completed.\n');
fprintf('Total time since computation started (including output): %.2f minutes (%.2f hours)\n', ...
    elapsedTotalSeconds / 60, elapsedTotalSeconds / 3600);

%% 4. LOCAL FUNCTIONS

end  % End of SKLGB_pilot_LSTMX main function

function [pred36, valPred36, nrmse3, tda3, bestP3, valBest3, ...
    rmse3, testRange3, success3, failureCode3, failureID3] = ...
    runLSTMXOnePixelAllRounds(X, Y, p_list, numRepeat, testYears, horizon, ...
    lstmParam, randomSeed, pointID)

    X = double(X);
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

    if ~isequal(size(X), [252, lstmParam.numFeatures]) || ...
            numel(Y) ~= 252 || any(~isfinite(X(:))) || any(~isfinite(Y))
        error('LSTMX input must be a finite 252-by-%d series.', ...
            lstmParam.numFeatures);
    end
    if any(X(:, 1) ~= Y)
        error('The first input column must equal the dGWS target series.');
    end

    for r = 1:nRound
        testYear = testYears(r);
        trainEndIdx = (testYear - 2005) * 12;
        valIdx = (trainEndIdx + 1):(trainEndIdx + horizon);
        testIdx = (trainEndIdx + horizon + 1): ...
            (trainEndIdx + 2 * horizon);
        predIdx = ((r - 1) * horizon + 1):(r * horizon);

        X_train = X(1:trainEndIdx, :);
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
            if size(X_train, 1) <= p + horizon - 1
                continue;
            end

            % Validation targets control early stopping and history selection.
            % No validation-year samples enter gradient-based training.
            X_valBlock = X(trainEndIdx-p+1:trainEndIdx+horizon, :);
            Y_valBlock = Y(trainEndIdx-p+1:trainEndIdx+horizon);
            initWindowVal = X_train(end-p+1:end, :);
            valScores = NaN(numRepeat, 1);
            candidatePacks = cell(numRepeat, 1);
            candidateValPred = NaN(horizon, numRepeat);

            for rep = 1:numRepeat
                try
                    seedVal = makeSeed(randomSeed, pointID, r, k, rep);
                    pack = fitLSTMXSeq2VectorFromSeries( ...
                        X_train, Y_train, p, lstmParam, seedVal, ...
                        X_valBlock, Y_valBlock);
                    predVal = predictLSTMXSeq2Vector(pack, initWindowVal);
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

            okRep = isfinite(valScores);
            if ~any(okRep)
                continue;
            end
            candidateScore = mean(valScores(okRep));
            if isfinite(candidateScore) && candidateScore < bestValScore
                bestValScore = candidateScore;
                bestP = p;
                selectedPacks = candidatePacks;
                selectedValPred = mean(candidateValPred(:, okRep), 2);
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

        % Testing updates the observed multivariate history, not the network.
        originIdx = trainEndIdx + horizon;
        initWindowTest = X(originIdx-bestP+1:originIdx, :);
        testPredRepeats = NaN(horizon, numRepeat);
        firstTestError = '';
        for rep = 1:numRepeat
            pack = selectedPacks{rep};
            if isempty(pack)
                continue;
            end
            try
                testPredRepeats(:, rep) = predictLSTMXSeq2Vector( ...
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

function lstmPack = fitLSTMXSeq2VectorFromSeries( ...
    X_series, Y_series, p, lstmParam, seedVal, X_valBlock, Y_valBlock)

    X_series = double(X_series);
    Y_series = double(Y_series(:));
    numPredict = lstmParam.numPredict;
    numFeatures = lstmParam.numFeatures;
    if size(X_series, 2) ~= numFeatures || ...
            size(X_series, 1) ~= numel(Y_series) || ...
            size(X_series, 1) <= p + numPredict - 1
        error('Invalid LSTMX training data or history length.');
    end

    [XTrainSeq, OutputTrain] = buildMultiSeq2VectorSamples( ...
        X_series, Y_series, p, numPredict);
    [XTrainNorm, mu_X, sigma_X] = standardizeTrainXSeq(XTrainSeq);
    [OutputTrainNorm, mu_Y, sigma_Y] = standardizeTrainY(OutputTrain);
    nTrain = size(XTrainNorm, 3);
    InputTrainCell = cell(nTrain, 1);
    for j = 1:nTrain
        InputTrainCell{j} = XTrainNorm(:, :, j);
    end

    layers = [
        sequenceInputLayer(numFeatures, 'Name', 'input')
        lstmLayer(lstmParam.numHiddenUnits, 'OutputMode', 'last', 'Name', 'lstm')
        fullyConnectedLayer(numPredict, 'Name', 'fc')
        regressionLayer('Name', 'regression')];
    validationFrequency = max(1, floor(nTrain / lstmParam.miniBatchSize));

    seedSafe = mod(seedVal, 2^32 - 1);
    if seedSafe == 0
        seedSafe = 1;
    end
    rng(seedSafe, 'twister');

    if ~isempty(X_valBlock)
        X_valBlock = double(X_valBlock);
        Y_valBlock = double(Y_valBlock(:));
        if ~isequal(size(X_valBlock), [p + numPredict, numFeatures]) || ...
                numel(Y_valBlock) ~= p + numPredict
            error('Validation block size must be [p+numPredict, numFeatures].');
        end
        InputVal = X_valBlock(1:p, :)';
        OutputVal = Y_valBlock(p+1:end)';
        InputValNorm = (InputVal - mu_X) ./ sigma_X;
        OutputValNorm = (OutputVal - mu_Y) ./ sigma_Y;
        InputValCell = {InputValNorm};
        options = trainingOptions(lstmParam.optimizer, ...
            'ExecutionEnvironment', lstmParam.executionEnvironment, ...
            'MaxEpochs', lstmParam.maxEpochs, ...
            'MiniBatchSize', lstmParam.miniBatchSize, ...
            'InitialLearnRate', lstmParam.initialLearnRate, ...
            'ValidationData', {InputValCell, OutputValNorm}, ...
            'ValidationFrequency', validationFrequency, ...
            'ValidationPatience', lstmParam.validationPatience, ...
            'OutputNetwork', lstmParam.outputNetwork, ...
            'Verbose', 0, ...
            'Plots', 'none', ...
            'Shuffle', lstmParam.shuffle, ...
            'GradientThreshold', lstmParam.gradientThreshold);
    else
        options = trainingOptions(lstmParam.optimizer, ...
            'ExecutionEnvironment', lstmParam.executionEnvironment, ...
            'MaxEpochs', lstmParam.maxEpochs, ...
            'MiniBatchSize', lstmParam.miniBatchSize, ...
            'InitialLearnRate', lstmParam.initialLearnRate, ...
            'Verbose', 0, ...
            'Plots', 'none', ...
            'Shuffle', lstmParam.shuffle, ...
            'GradientThreshold', lstmParam.gradientThreshold);
    end

    Mdl = trainNetwork(InputTrainCell, OutputTrainNorm, layers, options);
    lstmPack = struct();
    lstmPack.Mdl = Mdl;
    lstmPack.p = p;
    lstmPack.numPredict = numPredict;
    lstmPack.numFeatures = numFeatures;
    lstmPack.mu_X = mu_X;
    lstmPack.sigma_X = sigma_X;
    lstmPack.mu_Y = mu_Y;
    lstmPack.sigma_Y = sigma_Y;
end

function [X_seq, Y_target] = buildMultiSeq2VectorSamples( ...
    X_series, Y_series, p, numPredict)
    X_series = double(X_series);
    Y_series = double(Y_series(:));
    nTime = size(X_series, 1);
    nFeatures = size(X_series, 2);
    nSample = nTime - p - numPredict + 1;
    if nSample <= 0
        error('No sequence-to-vector samples can be constructed.');
    end
    X_seq = zeros(nFeatures, p, nSample);
    Y_target = zeros(nSample, numPredict);
    for j = 1:nSample
        X_seq(:, :, j) = X_series(j:j+p-1, :)';
        Y_target(j, :) = Y_series(j+p:j+p+numPredict-1)';
    end
end

function Y_pred = predictLSTMXSeq2Vector(lstmPack, initialWindow)
    p = lstmPack.p;
    initialWindow = double(initialWindow);
    if ~isequal(size(initialWindow), [p, lstmPack.numFeatures])
        error('Invalid LSTMX forecast-origin input window.');
    end
    inputSeq = initialWindow';
    inputNorm = (inputSeq - lstmPack.mu_X) ./ lstmPack.sigma_X;
    Y_pred_norm = predict(lstmPack.Mdl, inputNorm);
    Y_pred_norm = double(Y_pred_norm(:)');
    Y_pred = Y_pred_norm .* lstmPack.sigma_Y + lstmPack.mu_Y;
    Y_pred = double(Y_pred(:));
    if numel(Y_pred) ~= lstmPack.numPredict || any(~isfinite(Y_pred))
        error('LSTMX prediction is nonfinite or has the wrong length.');
    end
end

function [X_norm, mu_X, sigma_X] = standardizeTrainXSeq(X_seq)
    X_seq = double(X_seq);
    nFeatures = size(X_seq, 1);
    mu_X = zeros(nFeatures, 1);
    sigma_X = ones(nFeatures, 1);
    for f = 1:nFeatures
        tmp = X_seq(f, :, :);
        vals = tmp(:);
        mu_X(f) = mean(vals);
        sigma_X(f) = std(vals, 0);
        if ~isfinite(sigma_X(f)) || sigma_X(f) == 0
            sigma_X(f) = 1;
        end
    end
    X_norm = (X_seq - reshape(mu_X, [], 1, 1)) ./ ...
        reshape(sigma_X, [], 1, 1);
end

function [Y_norm, mu_Y, sigma_Y] = standardizeTrainY(Y)
    Y = double(Y);
    mu_Y = mean(Y, 1);
    sigma_Y = std(Y, 0, 1);
    sigma_Y(~isfinite(sigma_Y) | sigma_Y == 0) = 1;
    Y_norm = (Y - mu_Y) ./ sigma_Y;
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

function signature = computeRunSignature(dGWS, prec, evap, runoff, dSM2m, ...
    surfT, index, lat, lon, ym, ...
    p_list, testYears, horizon, numRepeat, lstmParam, randomSeed, blockNum, version)
    md = java.security.MessageDigest.getInstance('SHA-256');
    updateDigest(md, 'LSTMX_global_forecast');
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
    updateDigest(md, double(testYears));
    updateDigest(md, double(horizon));
    updateDigest(md, double(numRepeat));
    updateDigest(md, double(lstmParam.numFeatures));
    updateDigest(md, double(lstmParam.numHiddenUnits));
    updateDigest(md, double(lstmParam.initialLearnRate));
    updateDigest(md, double(lstmParam.miniBatchSize));
    updateDigest(md, double(lstmParam.maxEpochs));
    updateDigest(md, double(lstmParam.validationPatience));
    updateDigest(md, lstmParam.outputNetwork);
    updateDigest(md, double(lstmParam.gradientThreshold));
    updateDigest(md, lstmParam.optimizer);
    updateDigest(md, lstmParam.executionEnvironment);
    updateDigest(md, lstmParam.shuffle);
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
        'LSTMX_Pred_valid', 'LSTMX_Val_Pred_valid', ...
        'LSTMX_NRMSE_valid', 'LSTMX_TDA_valid', ...
        'LSTMX_best_p_valid', 'LSTMX_Val_NRMSE_best_valid', ...
        'LSTMX_RMSE_valid', 'LSTMX_TestRange_valid', 'LSTMX_TestSuccess_valid', ...
        'LSTMX_FailureCode_valid', 'LSTMX_FailureID_valid', ...
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
    matrix3 = {'LSTMX_NRMSE_valid', 'LSTMX_TDA_valid', ...
        'LSTMX_best_p_valid', 'LSTMX_Val_NRMSE_best_valid', ...
        'LSTMX_RMSE_valid', 'LSTMX_TestRange_valid', 'LSTMX_TestSuccess_valid', ...
        'LSTMX_FailureCode_valid', 'LSTMX_FailureID_valid'};
    for k = 1:numel(matrix3)
        if ~isequal(size(C.(matrix3{k})), [nValid, nRound])
            error('Checkpoint contains an invalid %s array.', matrix3{k});
        end
    end
    if ~isequal(size(C.LSTMX_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.LSTMX_Val_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.blockDone), [totalBlocks, 1])
        error('Checkpoint dimensions do not match the current run.');
    end
end
