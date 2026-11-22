%% SVMX - Global multivariate groundwater forecasting
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Forecast monthly groundwater-storage change (dGWS) at each valid
%   1-degree grid cell using support vector regression (SVR) and historical
%   hydroclimatic predictors.
%
% Input
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   Monthly [grid, month] series for dGWS, precipitation, evaporation,
%   runoff, soil-moisture change and surface temperature, 2004-2024.
%
% Rolling evaluation
%   Test 2022: train 2004-2020, validate 2021.
%   Test 2023: train 2004-2021, validate 2022.
%   Test 2024: train 2004-2022, validate 2023.
%   Each test forecast covers 12 months from the preceding December.
%
% Model
%   Gaussian (RBF) support vector regression, fitrsvm.
%   dGWS history length p: [1, 6, 12, 24, 36, 48, 60] months.
%   Exogenous history length x: [1, 3, 6] months for each of five variables.
%   BoxConstraint: [1, 2, 5, 10, 20, 50, 100].
%   KernelScale = 'auto'; Epsilon = 0.1 in standardized target units.
%   For each candidate (p, x, BoxConstraint), twelve separate regressors
%   predict horizons h = 1,...,12 from the same forecast origin.
%   All feature and target normalization parameters use training data only.
%   Validation-year NRMSE selects p, x and BoxConstraint. The selected fitted
%   regressors are reused for test prediction without refitting. Observed
%   validation-year data update the history at the test forecast origin.
%   No validation-year targets are included in parameter estimation, and no
%   exogenous data from the forecast year are used as predictors.
%
% Evaluation
%   NRMSE = RMSE / (max(test dGWS) - min(test dGWS)).
%   TDA compares the signs of successive month-to-month dGWS differences,
%   with 11 comparisons per test year.
%
% Output
%   SVMX_Pred, SVMX_Val_Pred: [lat, lon, 36].
%   SVMX_NRMSE, SVMX_TDA, SVMX_best_p, SVMX_best_x,
%   SVMX_best_box, SVMX_Val_NRMSE_best: [lat, lon, 3].
%   SVMX_RMSE, SVMX_TestRange, SVMX_TestSuccess: [lat, lon, 3].
%   SVMX_FailureCode: [lat, lon, 3]; code definitions in Meta.
%   Valid-grid arrays, spatial coordinates, dates and metadata are saved.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
function SKLGB_ext_SVMX(cfg)
% Original forecasting logic adapted to the fixed 2000-cell master.


projectRoot = cfg.projectRoot;
dataFile = cfg.pilotMasterFile;
outFolder = fullfile(cfg.outputRoot, 'SVMX');

% Independent experiments with distinct random seeds and output files.
experimentName = 'original';  % 'original', 'uc1', or 'uc2'
switch lower(experimentName)
    case 'original'
        randomSeed = 20260605;
        outFileName = 'SVMX_extended_2019_2021.mat';
        checkpointFileName = 'SVMX_extended_2019_2021_checkpoint.mat';
    case 'uc1'
        randomSeed = 20260719;
        outFileName = 'SVMX_multivariate_global_200401_202412_uc1.mat';
        checkpointFileName = 'SVMX_multivariate_global_uc1_checkpoint.mat';
    case 'uc2'
        randomSeed = 20260721;
        outFileName = 'SVMX_multivariate_global_200401_202412_uc2.mat';
        checkpointFileName = 'SVMX_multivariate_global_uc2_checkpoint.mat';
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
box_list = [1, 2, 5, 10, 20, 50, 100];
inputVarNames = {'dGWS', 'Prec', 'Evap', 'Runoff', 'dSM2m', 'SurfT'};
exogVarNames = {'Prec', 'Evap', 'Runoff', 'dSM2m', 'SurfT'};
testYears = [2019, 2020, 2021];
forecastHorizon = 12;

svmParam = struct();
svmParam.KernelFunction = 'rbf';
svmParam.KernelScale = 'auto';
svmParam.Epsilon = 0.1;

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
if any(box_list <= 0) || any(~isfinite(box_list)) || ...
        numel(unique(box_list)) ~= numel(box_list)
    error('box_list must contain distinct positive finite values.');
end
if ~isscalar(randomSeed) || randomSeed < 1 || ...
        randomSeed ~= floor(randomSeed) || randomSeed > 2^32 - 1
    error('randomSeed must be a positive integer within the uint32 range.');
end
if ~isequal(testYears, [2019, 2020, 2021]) || forecastHorizon ~= 12
    error('This experiment requires three annual 12-month forecasts, 2019-2021.');
end
if svmParam.Epsilon < 0 || ~isfinite(svmParam.Epsilon)
    error('SVMX Epsilon must be nonnegative and finite.');
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
originalValidRows = double(D.originalValidRows(:));
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

nRound = numel(testYears);
nPredMonth = nRound * forecastHorizon;
predYM = makeYMList(testYears(1), 1, testYears(end), 12);
valYM = makeYMList(testYears(1)-1, 1, testYears(end)-1, 12);
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);

block_size = max(1, ceil(nValid / block_num));
total_blocks = ceil(nValid / block_size);
fprintf('SVMX experiment: %s, random seed %d.\n', experimentName, randomSeed);
fprintf('Grid cells: %d; monthly observations: %d; test rounds: %d.\n', ...
    nValid, nTime, nRound);
fprintf('History lengths: %s\n', mat2str(p_list));
fprintf('Exogenous history lengths: %s\n', mat2str(x_list));
fprintf('Box constraints: %s\n', mat2str(box_list));
fprintf('Kernel: %s; kernel scale: %s; epsilon: %.3g.\n', ...
    svmParam.KernelFunction, svmParam.KernelScale, svmParam.Epsilon);
fprintf('Blocks: %d; approximately %d grid cells per block.\n', ...
    total_blocks, block_size);

%% 3.2 Prepare checkpoint and result arrays
runSignature = computeRunSignature(dGWS_valid, Prec_valid, Evap_valid, ...
    Runoff_valid, dSM2m_valid, SurfT_valid, validLinearIndex, ...
    lat, lon, ym, p_list, x_list, box_list, testYears, forecastHorizon, ...
    svmParam, randomSeed, block_num, checkpointVersion);

if isfile(checkpointFile) && ~resumeFromCheckpoint
    delete(checkpointFile);
end

if resumeFromCheckpoint && isfile(checkpointFile)
    fprintf('Reading checkpoint: %s\n', checkpointFile);
    C = load(checkpointFile);
    validateCheckpoint(C, runSignature, nValid, nRound, nPredMonth, total_blocks);

    SVMX_Pred_valid = C.SVMX_Pred_valid;
    SVMX_Val_Pred_valid = C.SVMX_Val_Pred_valid;
    SVMX_NRMSE_valid = C.SVMX_NRMSE_valid;
    SVMX_TDA_valid = C.SVMX_TDA_valid;
    SVMX_best_p_valid = C.SVMX_best_p_valid;
    SVMX_best_x_valid = C.SVMX_best_x_valid;
    SVMX_best_box_valid = C.SVMX_best_box_valid;
    SVMX_Val_NRMSE_best_valid = C.SVMX_Val_NRMSE_best_valid;
    SVMX_RMSE_valid = C.SVMX_RMSE_valid;
    SVMX_TestRange_valid = C.SVMX_TestRange_valid;
    SVMX_TestSuccess_valid = C.SVMX_TestSuccess_valid;
    SVMX_FailureCode_valid = C.SVMX_FailureCode_valid;
    SVMX_FailureID_valid = C.SVMX_FailureID_valid;
    blockDone = C.blockDone;
    cumulativeComputeTimeSeconds = C.cumulativeComputeTimeSeconds;
    clear C
else
    SVMX_Pred_valid = NaN(nValid, nPredMonth);
    SVMX_Val_Pred_valid = NaN(nValid, nPredMonth);
    SVMX_NRMSE_valid = NaN(nValid, nRound);
    SVMX_TDA_valid = NaN(nValid, nRound);
    SVMX_best_p_valid = NaN(nValid, nRound);
    SVMX_best_x_valid = NaN(nValid, nRound);
    SVMX_best_box_valid = NaN(nValid, nRound);
    SVMX_Val_NRMSE_best_valid = NaN(nValid, nRound);
    SVMX_RMSE_valid = NaN(nValid, nRound);
    SVMX_TestRange_valid = NaN(nValid, nRound);
    SVMX_TestSuccess_valid = false(nValid, nRound);
    SVMX_FailureCode_valid = uint8(zeros(nValid, nRound));
    SVMX_FailureID_valid = cell(nValid, nRound);
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

%% 3.4 Run global SVMX forecasts
fprintf('\nStarting SVMX computation...\n');
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
    blockBestBox = NaN(block_len, nRound);
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
            [p36, v36, n3, t3, best3, x3, box3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runSVMXOnePixelAllRounds(X_pixel, ...
                p_list, x_list, box_list, testYears, forecastHorizon, ...
                svmParam, randomSeed, originalValidRows(gridRow));
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockBestX(i, :) = x3(:)';
            blockBestBox(i, :) = box3(:)';
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
            [p36, v36, n3, t3, best3, x3, box3, vn3, rmse3, range3, ...
                ok3, code3, id3] = runSVMXOnePixelAllRounds(X_pixel, ...
                p_list, x_list, box_list, testYears, forecastHorizon, ...
                svmParam, randomSeed, originalValidRows(gridRow));
            blockPred(i, :) = p36(:)';
            blockValPred(i, :) = v36(:)';
            blockNRMSE(i, :) = n3(:)';
            blockTDA(i, :) = t3(:)';
            blockBestP(i, :) = best3(:)';
            blockBestX(i, :) = x3(:)';
            blockBestBox(i, :) = box3(:)';
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

    SVMX_Pred_valid(currentRows, :) = blockPred;
    SVMX_Val_Pred_valid(currentRows, :) = blockValPred;
    SVMX_NRMSE_valid(currentRows, :) = blockNRMSE;
    SVMX_TDA_valid(currentRows, :) = blockTDA;
    SVMX_best_p_valid(currentRows, :) = blockBestP;
    SVMX_best_x_valid(currentRows, :) = blockBestX;
    SVMX_best_box_valid(currentRows, :) = blockBestBox;
    SVMX_Val_NRMSE_best_valid(currentRows, :) = blockValBest;
    SVMX_RMSE_valid(currentRows, :) = blockRMSE;
    SVMX_TestRange_valid(currentRows, :) = blockTestRange;
    SVMX_TestSuccess_valid(currentRows, :) = blockSuccess;
    SVMX_FailureCode_valid(currentRows, :) = blockFailureCode;
    SVMX_FailureID_valid(currentRows, :) = blockFailureID;
    blockDone(block) = true;

    checkpointInfo = struct('version', checkpointVersion, ...
        'runSignature', runSignature);
    checkpointTemp = fullfile(outFolder, [checkpointFileName, '.tmp.mat']);
    save(checkpointTemp, ...
        'SVMX_Pred_valid', 'SVMX_Val_Pred_valid', ...
        'SVMX_NRMSE_valid', 'SVMX_TDA_valid', ...
        'SVMX_best_p_valid', 'SVMX_best_x_valid', 'SVMX_best_box_valid', ...
        'SVMX_Val_NRMSE_best_valid', 'SVMX_RMSE_valid', 'SVMX_TestRange_valid', ...
        'SVMX_TestSuccess_valid', 'SVMX_FailureCode_valid', 'SVMX_FailureID_valid', ...
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
fprintf('\nSVMX computation finished.\n');
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
SVMX_Pred = mapValidToGrid3D(SVMX_Pred_valid, validLinearIndex, nLat, nLon);
SVMX_Val_Pred = mapValidToGrid3D(SVMX_Val_Pred_valid, validLinearIndex, nLat, nLon);
SVMX_NRMSE = mapValidToGrid3D(SVMX_NRMSE_valid, validLinearIndex, nLat, nLon);
SVMX_TDA = mapValidToGrid3D(SVMX_TDA_valid, validLinearIndex, nLat, nLon);
SVMX_best_p = mapValidToGrid3D(SVMX_best_p_valid, validLinearIndex, nLat, nLon);
SVMX_best_x = mapValidToGrid3D(SVMX_best_x_valid, validLinearIndex, nLat, nLon);
SVMX_best_box = mapValidToGrid3D(SVMX_best_box_valid, validLinearIndex, nLat, nLon);
SVMX_Val_NRMSE_best = mapValidToGrid3D(SVMX_Val_NRMSE_best_valid, ...
    validLinearIndex, nLat, nLon);
SVMX_RMSE = mapValidToGrid3D(SVMX_RMSE_valid, validLinearIndex, nLat, nLon);
SVMX_TestRange = mapValidToGrid3D(SVMX_TestRange_valid, validLinearIndex, nLat, nLon);
SVMX_TestSuccess = mapValidToGrid3D(double(SVMX_TestSuccess_valid), ...
    validLinearIndex, nLat, nLon);
SVMX_FailureCode = mapValidToGrid3D(double(SVMX_FailureCode_valid), ...
    validLinearIndex, nLat, nLon);

SVMX_NRMSE_mean = mean(SVMX_NRMSE, 3, 'omitnan');
SVMX_TDA_mean = mean(SVMX_TDA, 3, 'omitnan');
SVMX_Val_NRMSE_best_mean = mean(SVMX_Val_NRMSE_best, 3, 'omitnan');

CompleteThreeYears_valid = all(SVMX_TestSuccess_valid, 2);
SVMX_CompleteThreeYears = false(nLat, nLon);
SVMX_CompleteThreeYears(validLinearIndex) = CompleteThreeYears_valid;
SVMX_NRMSE_mean_complete = SVMX_NRMSE_mean;
SVMX_NRMSE_mean_complete(~SVMX_CompleteThreeYears) = NaN;

fprintf('\nSVMX global results:\n');
fprintf('Successful grid-year forecasts: %d / %d.\n', ...
    nnz(SVMX_TestSuccess_valid), numel(SVMX_TestSuccess_valid));
fprintf('Complete three-year grid cells: %d / %d.\n', ...
    nnz(CompleteThreeYears_valid), nValid);
for r = 1:nRound
    fprintf('  %d success: %d / %d\n', testYears(r), ...
        nnz(SVMX_TestSuccess_valid(:, r)), nValid);
    for code = 1:4
        nFailure = nnz(SVMX_FailureCode_valid(:, r) == code);
        if nFailure > 0
            fprintf('    Failure code %d: %d\n', code, nFailure);
        end
    end
end
fprintf('SVMX_NRMSE_mean (finite values): %.4f\n', ...
    mean(SVMX_NRMSE_mean(CommonValid), 'omitnan'));
fprintf('SVMX_TDA_mean (finite values): %.2f%%\n', ...
    mean(SVMX_TDA_mean(CommonValid), 'omitnan'));

%% 3.6 Save final output
Meta = struct();
Meta.model = 'SVMX';
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
Meta.box_list = box_list;
Meta.svmParam = svmParam;
Meta.randomSeed = randomSeed;
Meta.standardization = ['Predictor columns and each horizon-specific target ' ...
    'are standardized using their respective training samples only'];
Meta.forecasting = 'Twelve direct SVR models, one for each lead month';
Meta.noFutureExogenousVariables = ['Only lagged exogenous variables that ' ...
    'are observed before each forecast origin are used'];
Meta.finalPrediction = ['Validation-selected fitted SVR models are reused ' ...
    'for testing without refitting; observed validation-year multivariate ' ...
    'history updates only predictors at the test forecast origin'];
Meta.randomization = ['KernelScale=auto uses randomized subsampling; ' ...
    'seeds depend on experiment, grid row, test round, p, x, box choice ' ...
    'and forecast horizon'];
Meta.NRMSE = 'Annual RMSE divided by observed test-year dGWS range';
Meta.TDA = ['Percentage of matching signs of successive month-to-month ' ...
    'differences in dGWS within a test year (11 comparisons)'];
Meta.failureCodes = ['0 = success; 1 = no finite validation candidate; ' ...
    '2 = test forecast error; 3 = invalid test predictions; ' ...
    '4 = nonfinite test NRMSE'];
Meta.SVMX_FailureID = 'MATLAB exception identifiers for failed fits or forecasts';
Meta.meanMetricRule = 'Three-year mean computed with omitnan';
Meta.completeMetricRule = ['SVMX_NRMSE_mean_complete requires successful ' ...
    'forecasts in all three test years'];
Meta.nValid = nValid;
Meta.block_num = block_num;
Meta.numWorkers = numWorkers;
Meta.usedParallel = useParallel;
Meta.runSignature = runSignature;
Meta.computeTimeMinutes = cumulativeComputeTimeSeconds / 60;
Meta.computeTimeThisRunMinutes = elapsedComputeThisRun / 60;
Meta.pilotOriginalValidRows = originalValidRows;
Meta.pilotMode = 'extended-2019-2021-original-hyperparameters';
Meta.scriptName = 'SKLGB_ext_SVMX.m';

finalTemp = fullfile(outFolder, [outFileName, '.tmp.mat']);
fprintf('\nSaving final MAT-file: %s\n', outFile);
save(finalTemp, ...
    'SVMX_Pred', 'SVMX_Val_Pred', ...
    'SVMX_NRMSE', 'SVMX_NRMSE_mean', 'SVMX_NRMSE_mean_complete', ...
    'SVMX_TDA', 'SVMX_TDA_mean', 'SVMX_best_p', 'SVMX_best_x', 'SVMX_best_box', ...
    'SVMX_Val_NRMSE_best', 'SVMX_Val_NRMSE_best_mean', ...
    'SVMX_RMSE', 'SVMX_TestRange', 'SVMX_TestSuccess', ...
    'SVMX_FailureCode', 'SVMX_CompleteThreeYears', ...
    'SVMX_Pred_valid', 'SVMX_Val_Pred_valid', ...
    'SVMX_NRMSE_valid', 'SVMX_TDA_valid', 'SVMX_best_p_valid', ...
    'SVMX_best_x_valid', ...
    'SVMX_best_box_valid', 'SVMX_Val_NRMSE_best_valid', ...
    'SVMX_RMSE_valid', 'SVMX_TestRange_valid', 'SVMX_TestSuccess_valid', ...
    'SVMX_FailureCode_valid', 'SVMX_FailureID_valid', ...
    'p_list', 'x_list', 'box_list', 'testYears', 'predYM', 'valYM', ...
    'svmParam', 'randomSeed', 'experimentName', ...
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
    figure('Name', 'SVMX mean NRMSE');
    p = imagesc(lon, lat, SVMX_NRMSE_mean);
    set(p, 'AlphaData', double(isfinite(SVMX_NRMSE_mean)));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(jet);
    clim(plotNRMSELimits);
    cb = colorbar;
    ylabel(cb, 'NRMSE');
    title('SVMX mean NRMSE, test years 2019-2021');
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

elapsedTotalSeconds = toc(totalTic);
fprintf('\nSVMX processing completed.\n');
fprintf('Total time since computation started (including output): %.2f minutes (%.2f hours)\n', ...
    elapsedTotalSeconds / 60, elapsedTotalSeconds / 3600);

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

end % End SKLGB_ext_SVMX main function

function [pred36, valPred36, nrmse3, tda3, bestP3, bestX3, bestBox3, ...
    valBest3, rmse3, testRange3, success3, failureCode3, failureID3] = ...
    runSVMXOnePixelAllRounds(X_all, p_list, x_list, box_list, testYears, ...
    horizon, svmParam, randomSeed, pointID)

    X_all = double(X_all);
    Y = X_all(:, 1);
    nRound = numel(testYears);
    pred36 = NaN(nRound * horizon, 1);
    valPred36 = NaN(nRound * horizon, 1);
    nrmse3 = NaN(nRound, 1);
    tda3 = NaN(nRound, 1);
    bestP3 = NaN(nRound, 1);
    bestX3 = NaN(nRound, 1);
    bestBox3 = NaN(nRound, 1);
    valBest3 = NaN(nRound, 1);
    rmse3 = NaN(nRound, 1);
    testRange3 = NaN(nRound, 1);
    success3 = false(nRound, 1);
    failureCode3 = uint8(zeros(nRound, 1));
    failureID3 = cell(nRound, 1);

    if ~isequal(size(X_all), [252, 6]) || any(~isfinite(X_all(:)))
        error('SVMX input must have six finite 252-month variables.');
    end

    for r = 1:nRound
        testYear = testYears(r);
        trainEndIdx = (testYear - 2005) * 12;
        valIdx = (trainEndIdx + 1):(trainEndIdx + horizon);
        testIdx = (trainEndIdx + horizon + 1): ...
            (trainEndIdx + 2 * horizon);
        predIdx = ((r - 1) * horizon + 1):(r * horizon);

        X_train = X_all(1:trainEndIdx, :);
        X_train_val = X_all(1:trainEndIdx + horizon, :);
        Y_train = Y(1:trainEndIdx);
        Y_val = Y(valIdx);
        Y_test = Y(testIdx);
        bestValScore = Inf;
        bestLinearIdx = Inf;
        selectedPack = [];
        valPredBest = [];
        firstValidationError = '';

        for ip = 1:numel(p_list)
            p = p_list(ip);
            for ix = 1:numel(x_list)
                x = x_list(ix);
                maxLag = max(p, x);
                if numel(Y_train) <= maxLag + horizon - 1
                    continue;
                end
                for ib = 1:numel(box_list)
                    boxVal = box_list(ib);
                    try
                        seedVal = makeSeed(randomSeed, pointID, r, ip, ix, ib);
                        svmPack = fitSVMXHorizonSpecificFromSeries( ...
                            X_train, Y_train, p, x, horizon, ...
                            svmParam, boxVal, seedVal);
                        initWindowVal = X_train(end-maxLag+1:end, :);
                        valPred = predictSVMXHorizonSpecific( ...
                            svmPack, initWindowVal);
                        score = calcNRMSE(Y_val, valPred);
                        candidateIdx = sub2ind( ...
                            [numel(p_list), numel(x_list), numel(box_list)], ...
                            ip, ix, ib);
                        if isfinite(score) && (score < bestValScore || ...
                                (score == bestValScore && ...
                                candidateIdx < bestLinearIdx))
                            bestValScore = score;
                            bestLinearIdx = candidateIdx;
                            selectedPack = svmPack;
                            valPredBest = valPred;
                        end
                    catch ME
                        if isempty(firstValidationError)
                            firstValidationError = ME.identifier;
                        end
                    end
                end
            end
        end

        if ~isfinite(bestValScore)
            failureCode3(r) = 1;
            failureID3{r} = firstValidationError;
            continue;
        end
        [bestIp, bestIx, bestIb] = ind2sub( ...
            [numel(p_list), numel(x_list), numel(box_list)], bestLinearIdx);
        bestP = p_list(bestIp);
        bestX = x_list(bestIx);
        bestBox = box_list(bestIb);
        bestP3(r) = bestP;
        bestX3(r) = bestX;
        bestBox3(r) = bestBox;
        valBest3(r) = bestValScore;
        valPred36(predIdx) = valPredBest(:);

        % Forecast from the new origin using the validation-selected models.
        % No fitted parameters or normalization statistics are recomputed.
        try
            maxLag = max(bestP, bestX);
            initWindowTest = X_train_val(end-maxLag+1:end, :);
            Y_test_pred = predictSVMXHorizonSpecific( ...
                selectedPack, initWindowTest);
        catch ME
            failureCode3(r) = 2;
            failureID3{r} = ME.identifier;
            continue;
        end

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

function svmPack = fitSVMXHorizonSpecificFromSeries( ...
    X_all, Y_series, p, x, horizon, svmParam, boxVal, baseSeed)

    modelCell = cell(horizon, 1);
    muXCell = cell(horizon, 1);
    sigmaXCell = cell(horizon, 1);
    muY = NaN(horizon, 1);
    sigmaY = NaN(horizon, 1);

    for h = 1:horizon
        [X, Y_target] = buildSVMXLaggedSamplesDirect( ...
            X_all, Y_series, p, x, h);
        [X_norm, mu_X, sigma_X] = standardizeTrainX(X);
        [Y_norm, mu_y, sigma_y] = standardizeTrainY(Y_target);
        rng(baseSeed + h, 'twister');
        mdl = fitrsvm(X_norm, Y_norm, ...
            'KernelFunction', svmParam.KernelFunction, ...
            'BoxConstraint', boxVal, ...
            'KernelScale', svmParam.KernelScale, ...
            'Epsilon', svmParam.Epsilon, ...
            'Standardize', false, ...
            'Verbose', 0);
        modelCell{h} = mdl;
        muXCell{h} = mu_X;
        sigmaXCell{h} = sigma_X;
        muY(h) = mu_y;
        sigmaY(h) = sigma_y;
    end

    svmPack = struct();
    svmPack.modelCell = modelCell;
    svmPack.p = p;
    svmPack.x = x;
    svmPack.horizon = horizon;
    svmPack.muXCell = muXCell;
    svmPack.sigmaXCell = sigmaXCell;
    svmPack.muY = muY;
    svmPack.sigmaY = sigmaY;
end

function [X, Y_target] = buildSVMXLaggedSamplesDirect( ...
    X_all, Y_series, p, x, h)

    nTime = size(X_all, 1);
    nExogVars = size(X_all, 2) - 1;
    maxLag = max(p, x);
    nSample = nTime - maxLag - h + 1;
    if nSample <= 0
        error('No direct SVMX training samples can be constructed.');
    end
    X = zeros(nSample, p + nExogVars * x);
    Y_target = zeros(nSample, 1);
    for j = 1:nSample
        anchorIdx = j + maxLag - 1;
        X(j, :) = makeSVMXOneFeatureRow( ...
            X_all(anchorIdx-maxLag+1:anchorIdx, :), p, x);
        Y_target(j) = Y_series(anchorIdx + h);
    end
end

function Y_pred = predictSVMXHorizonSpecific(svmPack, initialWindow)
    p = svmPack.p;
    x = svmPack.x;
    horizon = svmPack.horizon;
    maxLag = max(p, x);
    if size(initialWindow, 1) ~= maxLag || size(initialWindow, 2) ~= 6
        error('SVMX history window must have [max(p,x), 6] elements.');
    end
    X_new = makeSVMXOneFeatureRow(initialWindow, p, x);
    Y_pred = NaN(horizon, 1);
    for h = 1:horizon
        X_norm = (X_new - svmPack.muXCell{h}) ./ svmPack.sigmaXCell{h};
        y_norm = predict(svmPack.modelCell{h}, X_norm);
        Y_pred(h) = double(y_norm(1)) * svmPack.sigmaY(h) + svmPack.muY(h);
    end
end

function row = makeSVMXOneFeatureRow(initialWindow, p, x)
    nInputVars = size(initialWindow, 2);
    maxLag = max(p, x);
    if size(initialWindow, 1) ~= maxLag
        error('SVMX lag window has an unexpected number of time steps.');
    end
    row = zeros(1, p + (nInputVars - 1) * x);
    row(1:p) = initialWindow(end-p+1:end, 1)';
    pos = p + 1;
    for v = 2:nInputVars
        row(pos:pos+x-1) = initialWindow(end-x+1:end, v)';
        pos = pos + x;
    end
end

function [X_norm, mu_X, sigma_X] = standardizeTrainX(X)
    X = double(X);
    mu_X = mean(X, 1);
    sigma_X = std(X, 0, 1);
    sigma_X(~isfinite(sigma_X) | sigma_X == 0) = 1;
    X_norm = (X - mu_X) ./ sigma_X;
end

function [Y_norm, mu_Y, sigma_Y] = standardizeTrainY(Y)
    Y = double(Y(:));
    mu_Y = mean(Y);
    sigma_Y = std(Y, 0);
    if ~isfinite(sigma_Y) || sigma_Y == 0
        sigma_Y = 1;
    end
    Y_norm = (Y - mu_Y) ./ sigma_Y;
end

function seedVal = makeSeed(randomSeed, pointID, roundID, pID, xID, boxID)
    seedRaw = randomSeed + pointID * 100 + roundID * 10000 + ...
        pID * 1000 + xID * 100 + boxID * 10;
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
    surfT, index, lat, lon, ym, p_list, x_list, box_list, testYears, ...
    horizon, svmParam, randomSeed, blockNum, version)
    md = java.security.MessageDigest.getInstance('SHA-256');
    updateDigest(md, 'SVMX_global_forecast');
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
    updateDigest(md, double(box_list));
    updateDigest(md, double(testYears));
    updateDigest(md, double(horizon));
    updateDigest(md, svmParam.KernelFunction);
    updateDigest(md, svmParam.KernelScale);
    updateDigest(md, double(svmParam.Epsilon));
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
        'SVMX_Pred_valid', 'SVMX_Val_Pred_valid', ...
        'SVMX_NRMSE_valid', 'SVMX_TDA_valid', ...
        'SVMX_best_p_valid', 'SVMX_best_x_valid', 'SVMX_best_box_valid', ...
        'SVMX_Val_NRMSE_best_valid', 'SVMX_RMSE_valid', ...
        'SVMX_TestRange_valid', 'SVMX_TestSuccess_valid', ...
        'SVMX_FailureCode_valid', 'SVMX_FailureID_valid', ...
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
    matrix3 = {'SVMX_NRMSE_valid', 'SVMX_TDA_valid', ...
        'SVMX_best_p_valid', 'SVMX_best_x_valid', 'SVMX_best_box_valid', ...
        'SVMX_Val_NRMSE_best_valid', 'SVMX_RMSE_valid', ...
        'SVMX_TestRange_valid', 'SVMX_TestSuccess_valid', ...
        'SVMX_FailureCode_valid', 'SVMX_FailureID_valid'};
    for k = 1:numel(matrix3)
        if ~isequal(size(C.(matrix3{k})), [nValid, nRound])
            error('Checkpoint contains an invalid %s array.', matrix3{k});
        end
    end
    if ~isequal(size(C.SVMX_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.SVMX_Val_Pred_valid), [nValid, nPredMonth]) || ...
            ~isequal(size(C.blockDone), [totalBlocks, 1])
        error('Checkpoint dimensions do not match the current run.');
    end
end
