%% s5_NRMSE_sensitivity1.m
% ========================================================================
% Purpose:
%   Stage 1 of the NRMSE-normalization sensitivity analysis for manuscript 1.
%   Re-evaluate the EXISTING original-run predictions of nine models under
%   three normalization schemes, without retraining any forecasting model.
%
%   M0 = annual test RMSE / annual test observation range (current method)
%   A  = annual test RMSE / annual test observation IQR
%   B  = annual test RMSE / FIXED 2004-2020 observation IQR
%
%   For each grid cell and model, calculate each of the 2022, 2023, and
%   2024 annual scores separately and then take their arithmetic mean.
%   The denominator is common to the nine models at a given grid-year.
%
% Important scope:
%   - ORIGINAL model run only (not UC1 or UC2).
%   - No new model training, hyperparameter fitting, or forecasting.
%   - No winner selection or hierarchical attribution in this script;
%     they belong to s5_NRMSE_sensitivity2/3.
%   - Raw NRMSE numbers across M0/A/B are not on an identical scale.
%     Their model rankings and subsequent classifications can be compared.
%
% Reproducibility safeguards:
%   1. Test observations matched explicitly by ym = 202201 ... 202412.
%   2. Prediction layers interpreted as chronological 202201 ... 202412,
%      as in the existing s4_NRMSE_selection.m code.
%   3. Original common geographic mask and SAME complete-case pixel set
%      for the M0/A/B comparisons; no imputation or denominator clipping.
%   4. M0 is audited against saved MODEL_NRMSE and MODEL_NRMSE_mean maps.
%   5. IQR uses MATLAB's prctile(...,25/75,3) consistently for A and B.
%   6. Uses UNWEIGHTED global means, consistent with the manuscript.
%
% Inputs (revised master database):
%   03_master_database/GLDAS_common_valid_mask_1deg_200401_202412.mat
%     CommonValid; lat; lon
%   02_spatial_temporal_processing/
%     GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat
%     dGWS_mm; ym
%   Nine original-run model MAT files, each with MODEL_Pred (150x360x36)
%   and preferably MODEL_NRMSE, MODEL_NRMSE_mean for M0 audit.
%
% Output location:
%   <projectRoot>/s3_darft/s1_analysis/s5_NRMSE_sensitivity1_results/
%
% Outputs:
%   s5_NRMSE_sensitivity1_results.mat   (for stages 2 and 3)
%   s5_NRMSE_sensitivity1_tables.xlsx
%   s5_global_model_comparison.csv
%   s5_annual_model_comparison.csv
%   s5_M0_consistency_audit.csv
%   s5_denominator_diagnostics.csv
%   s5_mask_diagnostics.csv
%   s5_NRMSE_sensitivity1_log.txt
% ========================================================================

clc;
clearvars;

%% ===================== 1. Centralized settings =====================
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
expectedValidPixelCount = 17284;     % informative warning; not a mask rule
m0AuditTolerance = 1e-6;           % absolute NRMSE difference (audit)
smallDenominatorMM = 1e-6;         % DIAGNOSTIC only; no positive IQR clipped

historyStartYear = 2004;
historyEndYear   = 2020;           % fixed B reference; never expanded

testYears = [2022, 2023, 2024];
nYear = numel(testYears);

modelNames = { ...
    'AR'; 'ARIMA'; 'SARIMA'; ...
    'SVM'; 'RF'; 'LSTM'; ...
    'SVMX'; 'RFX'; 'LSTMX'};
nModel = numel(modelNames);

schemeNames = {'M0_test_range'; 'A_test_IQR'; 'B_fixed_history_IQR'};
nScheme = numel(schemeNames);

dataRoot = fullfile(projectRoot, 's1_master_database');
spatialTemporalDir = fullfile(dataRoot, '02_spatial_temporal_processing');
linearDir = fullfile(dataRoot, '04_linear_models');
nonlinearDir = fullfile(dataRoot, '05_nonlinear_models');
nonlinearXDir = fullfile(dataRoot, '06_nonlinear_x_models');

maskFile = fullfile(dataRoot, '03_master_database', ...
    'GLDAS_common_valid_mask_1deg_200401_202412.mat');
observationFile = fullfile(spatialTemporalDir, ...
    'GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat');

modelFiles = { ...
    fullfile(linearDir, 'AR_univariate_global_200401_202412.mat'); ...
    fullfile(linearDir, 'ARIMA_univariate_global_200401_202412_d0.mat'); ...
    fullfile(linearDir, 'SARIMA_univariate_global_200401_202412.mat'); ...
    fullfile(nonlinearDir, 'SVM_univariate_global_200401_202412.mat'); ...
    fullfile(nonlinearDir, 'RF_univariate_global_200401_202412.mat'); ...
    fullfile(nonlinearDir, 'LSTM_univariate_global_200401_202412.mat'); ...
    fullfile(nonlinearXDir, 'SVMX_multivariate_global_200401_202412.mat'); ...
    fullfile(nonlinearXDir, 'RFX_multivariate_global_200401_202412.mat'); ...
    fullfile(nonlinearXDir, 'LSTMX_multivariate_global_200401_202412.mat')};

outDir = fullfile(projectRoot, 's3_darft', 's1_analysis', ...
    's5_NRMSE_sensitivity1_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logFile = fullfile(outDir, 's5_NRMSE_sensitivity1_log.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
startTic = tic;
fprintf('============================================================\n');
fprintf('S5 STAGE 1: NRMSE NORMALIZATION SENSITIVITY (ORIGINAL RUN)\n');
fprintf('============================================================\n');
fprintf('M0: 12-month test Range | A: 12-month test IQR\n');
fprintf('B : fixed %d-%d historical IQR (all three test years)\n', ...
    historyStartYear, historyEndYear);
fprintf('Project root: %s\n\n', projectRoot);

%% ===================== 2. Check input files =====================
requiredFiles = [{maskFile; observationFile}; modelFiles];
for k = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{k})
        error('Required file is missing:\n%s', requiredFiles{k});
    end
end

%% ===================== 3. Read mask and observations =====================
S_mask = load(maskFile, 'CommonValid', 'lat', 'lon');
assert(isfield(S_mask,'CommonValid') && isfield(S_mask,'lat') && ...
       isfield(S_mask,'lon'), ...
       'Common mask MAT file must contain CommonValid, lat, lon.');

CommonValid = logical(S_mask.CommonValid);
lat = double(S_mask.lat(:));
lon = double(S_mask.lon(:));
[nLat, nLon] = size(CommonValid);
assert(nLat == numel(lat) && nLon == numel(lon), ...
    'CommonValid shape does not match lat/lon vectors.');
clear S_mask;

S_obs = load(observationFile, 'dGWS_mm', 'ym');
assert(isfield(S_obs,'dGWS_mm') && isfield(S_obs,'ym'), ...
    'Observation MAT file must contain dGWS_mm and ym.');

fullYM = double(S_obs.ym(:));
assert(numel(fullYM) == size(S_obs.dGWS_mm,3), ...
    'The number of ym dates must match dGWS_mm time layers.');
assert(numel(unique(fullYM)) == numel(fullYM), ...
    'The observation ym vector contains duplicate dates.');
assert(size(S_obs.dGWS_mm,1) == nLat && ...
       size(S_obs.dGWS_mm,2) == nLon, ...
       'Observation and valid-mask spatial dimensions differ.');

% Explicit, chronological YYYYMM vectors.
historyYM = buildMonthlyYM(historyStartYear, historyEndYear);
testYM = buildMonthlyYM(testYears(1), testYears(end));

[historyFound, historyIndices] = ismember(historyYM, fullYM);
[testFound, testIndices] = ismember(testYM, fullYM);
if ~all(historyFound)
    error('Missing history observation date(s), e.g. %d.', ...
        historyYM(find(~historyFound,1)));
end
if ~all(testFound)
    error('Missing test observation date(s), e.g. %d.', ...
        testYM(find(~testFound,1)));
end
assert(numel(testYM) == nYear*12, 'Expected 12 months per test year.');

observedHistory = double(S_obs.dGWS_mm(:,:,historyIndices));
observedTest = double(S_obs.dGWS_mm(:,:,testIndices));
clear S_obs;

fprintf('Original CommonValid pixels          : %d\n', nnz(CommonValid));
fprintf('Fixed historical months              : %d\n', numel(historyYM));
fprintf('Test months                          : %d\n', numel(testYM));
fprintf('Test years                           : %s\n\n', mat2str(testYears));

%% ===================== 4. Calculate three denominators =====================
% M0/A: each test year has its own 12-month observation reference.
% B: identical 2004-2020 IQR for every test year at the same grid cell.
% Neither IQR scheme is assumed to be superior to the original Range.

denom_M0_annual = nan(nLat,nLon,nYear);
denom_A_annual  = nan(nLat,nLon,nYear);

for iYear = 1:nYear
    monthIdx = (iYear-1)*12 + (1:12);
    observedYear = observedTest(:,:,monthIdx);

    denom_M0_annual(:,:,iYear) = ...
        max(observedYear,[],3) - min(observedYear,[],3);
    denom_A_annual(:,:,iYear) = ...
        prctile(observedYear,75,3) - prctile(observedYear,25,3);
end

% IQR over ALL 204 months in the fixed 2004-2020 history, not the
% test-year-specific expanding training windows.
denom_B_fixed = ...
    prctile(observedHistory,75,3) - prctile(observedHistory,25,3);

observedTestComplete = all(isfinite(observedTest),3);
observedHistoryComplete = all(isfinite(observedHistory),3);
clear observedHistory;

validM0Denom = all(isfinite(denom_M0_annual) & denom_M0_annual>0,3);
validADenom = all(isfinite(denom_A_annual) & denom_A_annual>0,3);
validBDenom = isfinite(denom_B_fixed) & denom_B_fixed>0;

fprintf('Complete 2022-2024 observed months   : %d within CommonValid\n', ...
    nnz(CommonValid & observedTestComplete));
fprintf('Complete 2004-2020 historical months: %d within CommonValid\n', ...
    nnz(CommonValid & observedHistoryComplete));
fprintf('Positive M0 annual ranges (all 3 y) : %d within CommonValid\n', ...
    nnz(CommonValid & validM0Denom));
fprintf('Positive A annual IQRs (all 3 y)    : %d within CommonValid\n', ...
    nnz(CommonValid & validADenom));
fprintf('Positive fixed B IQR                : %d within CommonValid\n\n', ...
    nnz(CommonValid & validBDenom));

%% ===================== 5. Calculate raw RMSE and M0/A/B =====================
% Array order for downstream scripts:
%   RMSE_annual / NRMSE_annual_* : [lat, lon, testYear, model]
%   NRMSE_mean_*                : [lat, lon, model]

RMSE_annual     = nan(nLat,nLon,nYear,nModel);
NRMSE_annual_M0 = nan(nLat,nLon,nYear,nModel);
NRMSE_annual_A  = nan(nLat,nLon,nYear,nModel);
NRMSE_annual_B  = nan(nLat,nLon,nYear,nModel);

% Optional saved-M0 comparison, populated when original variables exist.
storedMean_M0 = nan(nLat,nLon,nModel);
storedAnnual_M0 = nan(nLat,nLon,nYear,nModel);
hasStoredMean = false(nModel,1);
hasStoredAnnual = false(nModel,1);

modelPredictionComplete = true(nLat,nLon);

for iModel = 1:nModel
    thisFile = modelFiles{iModel};
    thisModel = modelNames{iModel};
    predVar = [thisModel, '_Pred'];
    meanVar = [thisModel, '_NRMSE_mean'];
    annualVar = [thisModel, '_NRMSE'];

    presentVars = who('-file', thisFile);
    if ~ismember(predVar, presentVars)
        error('Prediction variable %s is absent in:\n%s',predVar,thisFile);
    end

    S_model = load(thisFile, predVar);
    pred = double(S_model.(predVar));
    clear S_model;
    if ~isequal(size(pred),[nLat,nLon,numel(testYM)])
        error('%s shape mismatch: expected %s; got %s.', ...
            predVar,mat2str([nLat,nLon,numel(testYM)]),mat2str(size(pred)));
    end

    thisPredictionComplete = all(isfinite(pred),3);
    modelPredictionComplete = ...
        modelPredictionComplete & thisPredictionComplete;

    for iYear = 1:nYear
        monthIdx = (iYear-1)*12 + (1:12);
        err = pred(:,:,monthIdx) - observedTest(:,:,monthIdx);
        rmseMap = sqrt(mean(err.^2,3));
        RMSE_annual(:,:,iYear,iModel) = rmseMap;

        m0 = rmseMap ./ denom_M0_annual(:,:,iYear);
        a = rmseMap ./ denom_A_annual(:,:,iYear);
        b = rmseMap ./ denom_B_fixed;

        m0(~isfinite(m0)) = NaN;
        a(~isfinite(a)) = NaN;
        b(~isfinite(b)) = NaN;
        NRMSE_annual_M0(:,:,iYear,iModel) = m0;
        NRMSE_annual_A(:,:,iYear,iModel) = a;
        NRMSE_annual_B(:,:,iYear,iModel) = b;
    end

    clear pred err rmseMap m0 a b;

    if ismember(meanVar,presentVars)
        S_old = load(thisFile,meanVar);
        old = double(S_old.(meanVar));
        if isequal(size(old),[nLat,nLon])
            storedMean_M0(:,:,iModel) = old;
            hasStoredMean(iModel) = true;
        else
            warning('%s has unexpected dimensions and will not be audited.',meanVar);
        end
        clear old S_old;
    else
        warning('%s absent in %s; mean-level M0 audit not possible.', ...
            meanVar,thisFile);
    end

    if ismember(annualVar,presentVars)
        S_old = load(thisFile,annualVar);
        old = double(S_old.(annualVar));
        if isequal(size(old),[nLat,nLon,nYear])
            storedAnnual_M0(:,:,:,iModel) = old;
            hasStoredAnnual(iModel) = true;
        else
            warning('%s has unexpected dimensions and will not be audited.',annualVar);
        end
        clear old S_old;
    else
        warning('%s absent in %s; annual M0 audit not possible.', ...
            annualVar,thisFile);
    end

    fprintf('Loaded %d/%d : %-7s | finite prediction grids in CommonValid: %d\n', ...
        iModel,nModel,thisModel,nnz(CommonValid & thisPredictionComplete));
end

%% ===================== 6. Common comparison domain =====================
% The M0-only reference domain is retained for an audit against the
% manuscript data. All THREE sensitivity schemes use comparisonMask.
M0ReferenceMask = CommonValid & observedTestComplete & ...
    modelPredictionComplete & validM0Denom;
comparisonMask = M0ReferenceMask & observedHistoryComplete & ...
    validADenom & validBDenom;

nM0Ref = nnz(M0ReferenceMask);
nCommon = nnz(comparisonMask);
if nCommon == 0
    error('No grid cells remain valid across the three normalization schemes.');
end

fprintf('\n================ MASK AUDIT ================\n');
fprintf('CommonValid mask                       : %d\n',nnz(CommonValid));
fprintf('All-model M0 reference domain         : %d\n',nM0Ref);
fprintf('M0/A/B shared comparison domain       : %d\n',nCommon);
fprintf('Dropped from M0 domain for A/B        : %d\n',nM0Ref-nCommon);
if nCommon ~= expectedValidPixelCount
    warning(['Shared comparison domain has %d grids, rather than the ', ...
             'manuscript expected count %d. Inspect mask diagnostics.'], ...
            nCommon,expectedValidPixelCount);
end

maskCheck = { ...
    'CommonValid'; ...
    'Complete_test_observations'; ...
    'Complete_fixed_history_observations'; ...
    'Complete_all_model_predictions'; ...
    'Positive_M0_denominators_all_years'; ...
    'Positive_A_denominators_all_years'; ...
    'Positive_B_fixed_denominator'; ...
    'M0_reference_domain'; ...
    'M0_A_B_common_comparison_domain'; ...
    'Removed_by_A_B_relative_to_M0'};
maskCount = [ ...
    nnz(CommonValid); ...
    nnz(CommonValid & observedTestComplete); ...
    nnz(CommonValid & observedHistoryComplete); ...
    nnz(CommonValid & modelPredictionComplete); ...
    nnz(CommonValid & validM0Denom); ...
    nnz(CommonValid & validADenom); ...
    nnz(CommonValid & validBDenom); ...
    nM0Ref; nCommon; nM0Ref-nCommon];
maskTable = table(maskCheck,maskCount, ...
    'VariableNames',{'Check','PixelCount'});

%% ===================== 7. Reproduce original M0 scores =====================
% Do this BEFORE applying the A/B common-domain restrictions. Annual and
% mean maps are compared with their saved original-run counterparts.
% A mismatch is flagged, not hidden by recomputing or changing inputs.

% Calculate the three-year MEAN OF ANNUAL scores (not a 36-month RMSE).
NRMSE_mean_M0 = nan(nLat,nLon,nModel);
NRMSE_mean_A  = nan(nLat,nLon,nModel);
NRMSE_mean_B  = nan(nLat,nLon,nModel);

for iModel = 1:nModel
    NRMSE_mean_M0(:,:,iModel) = mean(NRMSE_annual_M0(:,:,:,iModel),3);
    NRMSE_mean_A(:,:,iModel) = mean(NRMSE_annual_A(:,:,:,iModel),3);
    NRMSE_mean_B(:,:,iModel) = mean(NRMSE_annual_B(:,:,:,iModel),3);
end

oldGlobalMeanOriginalMask = nan(nModel,1);
newGlobalMeanM0Ref = nan(nModel,1);
oldGlobalMeanM0Ref = nan(nModel,1);
maxAbsMeanDiff = nan(nModel,1);
maxAbsAnnualDiff = nan(nModel,1);
nComparedPixels = zeros(nModel,1);
auditStatus = repmat({'NOT_AVAILABLE'},nModel,1);

fprintf('\n================ ORIGINAL M0 REPRODUCTION AUDIT ================\n');
for iModel = 1:nModel
    newMap = NRMSE_mean_M0(:,:,iModel);
    newGlobalMeanM0Ref(iModel) = mean(newMap(M0ReferenceMask));

    if hasStoredMean(iModel)
        oldMap = storedMean_M0(:,:,iModel);
        if all(isfinite(oldMap(CommonValid)))
            oldGlobalMeanOriginalMask(iModel) = mean(oldMap(CommonValid));
        end

        auditMask = M0ReferenceMask & isfinite(oldMap) & isfinite(newMap);
        nComparedPixels(iModel) = nnz(auditMask);
        if any(auditMask(:))
            oldGlobalMeanM0Ref(iModel) = mean(oldMap(auditMask));
            maxAbsMeanDiff(iModel) = max(abs(oldMap(auditMask)-newMap(auditMask)));
        end
    end

    if hasStoredAnnual(iModel)
        diffYearMax = nan(nYear,1);
        for iYear = 1:nYear
            oldAnn = storedAnnual_M0(:,:,iYear,iModel);
            newAnn = NRMSE_annual_M0(:,:,iYear,iModel);
            auditYearMask = M0ReferenceMask & isfinite(oldAnn) & isfinite(newAnn);
            if any(auditYearMask(:))
                diffYearMax(iYear) = ...
                    max(abs(oldAnn(auditYearMask)-newAnn(auditYearMask)));
            end
        end
        if all(isfinite(diffYearMax))
            maxAbsAnnualDiff(iModel) = max(diffYearMax);
        end
    end

    if hasStoredMean(iModel) && hasStoredAnnual(iModel)
        if nComparedPixels(iModel) == nM0Ref && ...
                isfinite(maxAbsMeanDiff(iModel)) && ...
                isfinite(maxAbsAnnualDiff(iModel)) && ...
                maxAbsMeanDiff(iModel) <= m0AuditTolerance && ...
                maxAbsAnnualDiff(iModel) <= m0AuditTolerance
            auditStatus{iModel} = 'PASS';
        else
            auditStatus{iModel} = 'WARN_MISMATCH';
        end
    end
    fprintf('%-7s | old=%.8f new(M0-ref)=%.8f | max pixel diff=%.3g', ...
        modelNames{iModel},oldGlobalMeanM0Ref(iModel), ...
        newGlobalMeanM0Ref(iModel),maxAbsMeanDiff(iModel));
    fprintf(' | max year diff=%.3g | %s\n', ...
        maxAbsAnnualDiff(iModel),auditStatus{iModel});
end

M0AuditTable = table(modelNames,hasStoredMean,hasStoredAnnual, ...
    nComparedPixels,oldGlobalMeanOriginalMask,oldGlobalMeanM0Ref, ...
    newGlobalMeanM0Ref,maxAbsMeanDiff,maxAbsAnnualDiff,auditStatus, ...
    'VariableNames',{'Model','HasStoredMean','HasStoredAnnual', ...
    'ComparedPixels','OldGlobalMeanOriginalMask','OldGlobalMeanM0Domain', ...
    'RecomputedGlobalMeanM0Domain','MaxAbsMeanMapDifference', ...
    'MaxAbsAnnualDifference','AuditStatus'});
if any(strcmp(auditStatus,'WARN_MISMATCH'))
    warning(['At least one model failed the M0 reproduction audit. ', ...
             'Do not interpret A/B sensitivity until this is resolved.']);
end

%% ===================== 8. Enforce identical domain M0/A/B =====================
for iModel = 1:nModel
    mean0 = NRMSE_mean_M0(:,:,iModel);
    meanA = NRMSE_mean_A(:,:,iModel);
    meanB = NRMSE_mean_B(:,:,iModel);
    mean0(~comparisonMask) = NaN;
    meanA(~comparisonMask) = NaN;
    meanB(~comparisonMask) = NaN;
    NRMSE_mean_M0(:,:,iModel) = mean0;
    NRMSE_mean_A(:,:,iModel) = meanA;
    NRMSE_mean_B(:,:,iModel) = meanB;

    for iYear = 1:nYear
        rmse = RMSE_annual(:,:,iYear,iModel);
        m0 = NRMSE_annual_M0(:,:,iYear,iModel);
        a  = NRMSE_annual_A(:,:,iYear,iModel);
        b  = NRMSE_annual_B(:,:,iYear,iModel);
        rmse(~comparisonMask) = NaN;
        m0(~comparisonMask) = NaN;
        a(~comparisonMask) = NaN;
        b(~comparisonMask) = NaN;
        RMSE_annual(:,:,iYear,iModel) = rmse;
        NRMSE_annual_M0(:,:,iYear,iModel) = m0;
        NRMSE_annual_A(:,:,iYear,iModel) = a;
        NRMSE_annual_B(:,:,iYear,iModel) = b;
    end
end

%% ===================== 9. Global and annual comparisons =====================
% globalMeanScores/model ranks always refer to SAME grid-cell set.
globalMeanScores = nan(nModel,nScheme);
globalMedianScores = nan(nModel,nScheme);
globalRanks = nan(nModel,nScheme);
globalMeanRMSE = nan(nModel,1);
annualGlobalMeanScores = nan(nYear,nModel,nScheme);
annualGlobalMeanRMSE = nan(nYear,nModel);
annualRanks = nan(nYear,nModel,nScheme);

for iModel = 1:nModel
    scoreMaps = cat(4, ...
        NRMSE_mean_M0(:,:,iModel), ...
        NRMSE_mean_A(:,:,iModel), ...
        NRMSE_mean_B(:,:,iModel));

    for iScheme = 1:nScheme
        map = scoreMaps(:,:,1,iScheme);
        validScores = map(comparisonMask);
        if ~all(isfinite(validScores))
            error('Nonfinite mean NRMSE within comparisonMask (%s, %s).', ...
                modelNames{iModel},schemeNames{iScheme});
        end
        globalMeanScores(iModel,iScheme) = mean(validScores);
        globalMedianScores(iModel,iScheme) = median(validScores);
    end

    rmseMeanMap = mean(RMSE_annual(:,:,:,iModel),3);
    globalMeanRMSE(iModel) = mean(rmseMeanMap(comparisonMask));

    for iYear = 1:nYear
        annualRMSEMap = RMSE_annual(:,:,iYear,iModel);
        annualGlobalMeanRMSE(iYear,iModel) = ...
            mean(annualRMSEMap(comparisonMask));
        annualScores = cat(4, ...
            NRMSE_annual_M0(:,:,iYear,iModel), ...
            NRMSE_annual_A(:,:,iYear,iModel), ...
            NRMSE_annual_B(:,:,iYear,iModel));
        for iScheme = 1:nScheme
            map = annualScores(:,:,1,iScheme);
            vals = map(comparisonMask);
            if ~all(isfinite(vals))
                error('Nonfinite annual NRMSE in common mask (%s, %d, %s).', ...
                    modelNames{iModel},testYears(iYear),schemeNames{iScheme});
            end
            annualGlobalMeanScores(iYear,iModel,iScheme) = mean(vals);
        end
    end
end

% Deterministic ranking: ties broken by the defined nine-model order.
for iScheme = 1:nScheme
    [~,order] = sortrows([globalMeanScores(:,iScheme),(1:nModel)'],[1,2]);
    globalRanks(order,iScheme) = (1:nModel)';
    for iYear = 1:nYear
        yearScores = reshape(annualGlobalMeanScores(iYear,:,iScheme),nModel,1);
        [~,yearOrder] = sortrows([yearScores,(1:nModel)'],[1,2]);
        annualRanks(iYear,yearOrder,iScheme) = 1:nModel;
    end
end

globalTable = table(modelNames,globalMeanRMSE, ...
    globalMeanScores(:,1),globalRanks(:,1), ...
    globalMeanScores(:,2),globalRanks(:,2), ...
    globalMeanScores(:,3),globalRanks(:,3), ...
    'VariableNames',{'Model','MeanAnnualRMSE_mm', ...
    'M0_GlobalMeanNRMSE','M0_Rank', ...
    'A_GlobalMeanNRMSE','A_Rank', ...
    'B_GlobalMeanNRMSE','B_Rank'});

annualModelCol = cell(nModel*nYear,1);
annualYearCol = zeros(nModel*nYear,1);
annualRMSECol = nan(nModel*nYear,1);
annualM0Col = nan(nModel*nYear,1);
annualACol = nan(nModel*nYear,1);
annualBCol = nan(nModel*nYear,1);
annualM0RankCol = nan(nModel*nYear,1);
annualARankCol = nan(nModel*nYear,1);
annualBRankCol = nan(nModel*nYear,1);

row = 0;
for iYear = 1:nYear
    for iModel = 1:nModel
        row = row+1;
        annualModelCol{row} = modelNames{iModel};
        annualYearCol(row) = testYears(iYear);
        annualRMSECol(row) = annualGlobalMeanRMSE(iYear,iModel);
        annualM0Col(row) = annualGlobalMeanScores(iYear,iModel,1);
        annualACol(row)  = annualGlobalMeanScores(iYear,iModel,2);
        annualBCol(row)  = annualGlobalMeanScores(iYear,iModel,3);
        annualM0RankCol(row) = annualRanks(iYear,iModel,1);
        annualARankCol(row)  = annualRanks(iYear,iModel,2);
        annualBRankCol(row)  = annualRanks(iYear,iModel,3);
    end
end

annualTable = table(annualYearCol,annualModelCol,annualRMSECol, ...
    annualM0Col,annualM0RankCol,annualACol,annualARankCol, ...
    annualBCol,annualBRankCol, ...
    'VariableNames',{'TestYear','Model','GlobalMeanRMSE_mm', ...
    'M0_GlobalMeanNRMSE','M0_Rank','A_GlobalMeanNRMSE','A_Rank', ...
    'B_GlobalMeanNRMSE','B_Rank'});

%% ===================== 10. Denominator diagnostics =====================
% Report the denominator distributions WITHIN the shared analysis domain.
% Also report how many positive denominators are unusually small in mm.
% These diagnostics are part of this run, not a separate experiment.

diagScheme = cell(2*nYear+1,1);
diagYear = nan(2*nYear+1,1);
diagCount = zeros(2*nYear+1,1);
diagSmallCount = zeros(2*nYear+1,1);
diagMin = nan(2*nYear+1,1);
diagP01 = nan(2*nYear+1,1);
diagP05 = nan(2*nYear+1,1);
diagP25 = nan(2*nYear+1,1);
diagMedian = nan(2*nYear+1,1);
diagP75 = nan(2*nYear+1,1);
diagP95 = nan(2*nYear+1,1);
diagP99 = nan(2*nYear+1,1);
diagMax = nan(2*nYear+1,1);

for row = 1:(2*nYear+1)
    if row <= nYear
        thisDenom = denom_M0_annual(:,:,row);
        diagScheme{row} = schemeNames{1};
        diagYear(row) = testYears(row);
    elseif row <= 2*nYear
        iYear = row-nYear;
        thisDenom = denom_A_annual(:,:,iYear);
        diagScheme{row} = schemeNames{2};
        diagYear(row) = testYears(iYear);
    else
        thisDenom = denom_B_fixed;
        diagScheme{row} = schemeNames{3};
        % NaN year means one fixed 2004-2020 historical reference.
    end
    values = thisDenom(comparisonMask);
    diagCount(row) = numel(values);
    diagSmallCount(row) = nnz(values <= smallDenominatorMM);
    diagMin(row) = min(values);
    qq = prctile(values,[1,5,25,50,75,95,99]);
    diagP01(row) = qq(1);
    diagP05(row) = qq(2);
    diagP25(row) = qq(3);
    diagMedian(row) = qq(4);
    diagP75(row) = qq(5);
    diagP95(row) = qq(6);
    diagP99(row) = qq(7);
    diagMax(row) = max(values);
end

denominatorTable = table(diagScheme,diagYear,diagCount, ...
    diagSmallCount,diagMin,diagP01,diagP05,diagP25, ...
    diagMedian,diagP75,diagP95,diagP99,diagMax, ...
    'VariableNames',{'Scheme','TestYear','PixelCount', ...
    'NumDenomLE_1eMinus6mm','Min_mm','P01_mm','P05_mm', ...
    'P25_mm','Median_mm','P75_mm','P95_mm','P99_mm','Max_mm'});

%% ===================== 11. Save outputs =====================
% Data arrays are not predictions. They can be passed into stage 2 without
% reloading any model MAT file.

matFile = fullfile(outDir,'s5_NRMSE_sensitivity1_results.mat');
save(matFile, ...
    'modelNames','modelFiles','schemeNames', ...
    'testYears','testYM','historyYM','lat','lon', ...
    'projectRoot','maskFile','observationFile', ...
    'CommonValid','M0ReferenceMask','comparisonMask', ...
    'denom_M0_annual','denom_A_annual','denom_B_fixed', ...
    'RMSE_annual', ...
    'NRMSE_annual_M0','NRMSE_annual_A','NRMSE_annual_B', ...
    'NRMSE_mean_M0','NRMSE_mean_A','NRMSE_mean_B', ...
    'globalMeanScores','globalMedianScores','globalRanks', ...
    'annualGlobalMeanScores','annualRanks', ...
    'globalTable','annualTable','M0AuditTable', ...
    'maskTable','denominatorTable', ...
    'm0AuditTolerance','smallDenominatorMM', ...
    '-v7.3');

writetable(globalTable, ...
    fullfile(outDir,'s5_global_model_comparison.csv'));
writetable(annualTable, ...
    fullfile(outDir,'s5_annual_model_comparison.csv'));
writetable(M0AuditTable, ...
    fullfile(outDir,'s5_M0_consistency_audit.csv'));
writetable(maskTable, ...
    fullfile(outDir,'s5_mask_diagnostics.csv'));
writetable(denominatorTable, ...
    fullfile(outDir,'s5_denominator_diagnostics.csv'));

excelFile = fullfile(outDir,'s5_NRMSE_sensitivity1_tables.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(globalTable,excelFile,'Sheet','GlobalModels');
writetable(annualTable,excelFile,'Sheet','AnnualModels');
writetable(M0AuditTable,excelFile,'Sheet','M0Audit');
writetable(maskTable,excelFile,'Sheet','MaskSummary');
writetable(denominatorTable,excelFile,'Sheet','Denominators');

fprintf('\n================ GLOBAL MODEL COMPARISON ================\n');
disp(globalTable);
fprintf('\n================ YEARLY PERFORMANCE ================\n');
disp(annualTable);
fprintf('\nTotal time: %.1f seconds\n',toc(startTic));
fprintf('MAT result: %s\n',matFile);
fprintf('Excel workbook: %s\n',excelFile);

if all(strcmp(auditStatus,'PASS'))
    fprintf('\nM0 REPRODUCTION CHECK: ALL NINE MODELS PASS.\n');
else
    fprintf(['\nWARNING: M0 audit is not PASS for all models. ', ...
        'Review s5_M0_consistency_audit.csv before interpreting A/B.\n']);
end
if nCommon ~= nM0Ref
    fprintf(['WARNING: M0/A/B share %d pixels, but the M0-only domain has %d. ', ...
        'Review s5_mask_diagnostics.csv.\n'],nCommon,nM0Ref);
end
fprintf('============================================================\n');

diary off;

%% ===================== Local helper =====================
function monthlyYM = buildMonthlyYM(firstYear,lastYear)
    if firstYear > lastYear
        error('buildMonthlyYM requires firstYear <= lastYear.');
    end
    nYearLocal = lastYear-firstYear+1;
    monthlyYM = zeros(nYearLocal*12,1);
    for iy = 1:nYearLocal
        y = firstYear+iy-1;
        monthlyYM((iy-1)*12+(1:12)) = y*100+(1:12)';
    end
end
