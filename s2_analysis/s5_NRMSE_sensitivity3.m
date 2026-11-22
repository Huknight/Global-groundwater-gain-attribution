%% s5_NRMSE_sensitivity3.m
% =========================================================================
% STAGE 3: Hierarchical gain-COMPONENT sensitivity to NRMSE normalization
%
% INPUT (from Stage 2; NO model retraining or re-forecasting):
%   s5_NRMSE_sensitivity2_results/s5_NRMSE_sensitivity2_results.mat
%
% INPUT (from Stage 1; denominator diagnostics only):
%   s5_NRMSE_sensitivity1_results/s5_NRMSE_sensitivity1_results.mat
%
% Three evaluated NRMSE definitions:
%   M0: test-year RMSE / test-year observed Range (original manuscript)
%   A : test-year RMSE / test-year observed IQR
%   B : test-year RMSE / fixed 2004-2020 observed IQR
%
% This code exactly follows the six-class, winner-conditioned component
% definitions in s2_gain_attribution3.m and Appendix D of the manuscript.
% It uses the Stage 2 group references (L, N, X), winner, and classes.
%
% MAIN statistics always use ALL 17,284 common valid grid cells. The Table 1
% shares are ratios of SUMMED component reductions to SUMMED total reductions
% within each category. They are NOT means of pixel-specific percentages.
%
% Separate, labeled diagnostics check whether a small number of cells can
% dominate those summed shares. The diagnostic exclusions do NOT change
% the main results, winner identities, category composition, or study mask.
% Each exclusion is independently defined within a scheme and category.
%
% OUTPUTS (new directory; originals and earlier stages are left untouched):
%   s5_NRMSE_sensitivity3_results.mat
%   s5_NRMSE_sensitivity3_tables.xlsx
%   s5_stage3_Table1_M0.csv / _A.csv / _B.csv
%   s5_stage3_Hydro_main_comparison.csv
%   s5_stage3_component_shares_all_categories.csv
%   s5_stage3_category_statistics.csv
%   s5_stage3_extreme_cell_diagnostics.csv
%   s5_stage3_extreme_cell_concentration.csv
%   s5_stage3_M0_audit.csv
%   Fig1_Hydro_Contribution_Sensitivity.png (optional)
%   s5_NRMSE_sensitivity3_log.txt
%
% NOTE: 71% (history: Linear + Nonlinear) and 29% (Hydroclimatic) apply
% ONLY WITHIN THE HYDROCLIMATIC-GAIN CATEGORY, not across all 17,284 cells.
% =========================================================================

clc;
clearvars;

%% =================== 1. Centralized settings ==========================
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
baseDir = fullfile(projectRoot, 's3_darft', 's1_analysis');

stage2File = fullfile(baseDir, 's5_NRMSE_sensitivity2_results', ...
    's5_NRMSE_sensitivity2_results.mat');
stage1File = fullfile(baseDir, 's5_NRMSE_sensitivity1_results', ...
    's5_NRMSE_sensitivity1_results.mat');

outDir = fullfile(baseDir, 's5_NRMSE_sensitivity3_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Diagnostics only. Each category is examined separately.
% The main decomposition always includes every valid member of its class.
extremeFraction = 0.01;           % Top/bottom 1% within each class.
expectedValidCells = 17284;       % Informational warning, not a filter.
makeFigure = true;

% Audit M0 against Table 1 in the current manuscript. The tolerances match
% the precision of printed values: one decimal place for percentages and
% four decimal places for mean component reductions.
auditShareTolPP = 0.051;
auditMeanTol = 0.000051;
expectedM0Counts = [1406;4343;4503;2083;3008;1941];
expectedM0Shares = [ ...
    NaN,  NaN,  NaN,  NaN; ...  % No gain: zero denominator
    100,  0,    0,    0; ...
    30.8, 69.2, 0,    0; ...
    26.3, 44.7, 29.0, 0; ...
    58.3, 0,    0,    41.7; ...
    24.2, 45.6, 0,    30.2];
expectedM0MeanComponents = [ ...
    0,      0,      0,      0; ...
    0.0310, 0,      0,      0; ...
    0.0134, 0.0301, 0,      0; ...
    0.0118, 0.0201, 0.0130, 0; ...
    0.0167, 0,      0,      0.0120; ...
    0.0108, 0.0203, 0,      0.0135];

modelNamesExpected = { ...
    'AR';'ARIMA';'SARIMA';'SVM';'RF';'LSTM';'SVMX';'RFX';'LSTMX'};
schemeNamesExpected = { ...
    'M0_test_range';'A_test_IQR';'B_fixed_history_IQR'};
classNamesExpected = { ...
    'No gain';'Linear gain';'Nonlinear gain';'Hydroclimatic gain'; ...
    'Mixed-L';'Mixed-N'};
componentNames = {'Linear','Nonlinear','Hydroclimatic','Mixed'};
shortSchemeNames = {'M0';'A';'B'};

nScheme = numel(schemeNamesExpected);
nClass = numel(classNamesExpected);
hydroClassIdx = 4;

logFile = fullfile(outDir, 's5_NRMSE_sensitivity3_log.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
startTic = tic;

fprintf('===============================================================\n');
fprintf('S5 STAGE 3: HIERARCHICAL GAIN COMPONENT SENSITIVITY\n');
fprintf('===============================================================\n');
fprintf('M0 = test-year Range | A = test-year IQR | B = fixed historical IQR\n');
fprintf('Stage 2 input: %s\n', stage2File);
fprintf('Stage 1 denominator input: %s\n\n', stage1File);

%% ======================= 2. Read Stage 2 ============================
if ~isfile(stage2File)
    error('Stage 2 MAT file not found:\n%s', stage2File);
end
if ~isfile(stage1File)
    error('Stage 1 MAT file not found:\n%s', stage1File);
end

requiredStage2 = { ...
    'R_M0','R_A','R_B', ...
    'AR_NRMSE_M0','AR_NRMSE_A','AR_NRMSE_B', ...
    'modelNames','schemeNames','categoryNames','comparisonMask', ...
    'categoryCounts','testYears','lat','lon'};
checkMATVariables(stage2File,requiredStage2);
S = load(stage2File,requiredStage2{:});

modelNames = S.modelNames(:);
schemeNames = S.schemeNames(:);
classNames = S.categoryNames(:);
testYears = double(S.testYears(:)');
lat = double(S.lat(:));
lon = double(S.lon(:));
comparisonMask = logical(S.comparisonMask);
originalCategoryCounts = double(S.categoryCounts);

assert(isequal(modelNames,modelNamesExpected), ...
    'Unexpected nine-model order in Stage 2.');
assert(isequal(schemeNames,schemeNamesExpected), ...
    'Unexpected M0/A/B scheme names or order in Stage 2.');
assert(isequal(classNames,classNamesExpected), ...
    'Unexpected six-category names or order in Stage 2.');
assert(isequal(testYears,[2022 2023 2024]), ...
    'Unexpected test-year sequence in Stage 2.');
assert(isequal(size(comparisonMask),[numel(lat),numel(lon)]), ...
    'Stage 2 comparison mask and lat/lon dimensions disagree.');
assert(isequal(size(originalCategoryCounts),[nClass,nScheme]), ...
    'Stage 2 category count matrix should be 6-by-3.');

nValid = nnz(comparisonMask);
if nValid ~= expectedValidCells
    warning('Shared comparison mask has %d cells, expected %d.', ...
        nValid,expectedValidCells);
end

R_input = {S.R_M0;S.R_A;S.R_B};
AR_input = {double(S.AR_NRMSE_M0); ...
            double(S.AR_NRMSE_A); ...
            double(S.AR_NRMSE_B)};
clear S;

fprintf('Shared valid grid cells: %d\n',nValid);
fprintf('Test years: %s\n\n',mat2str(testYears));

%% =============== 3. Read denominator diagnostic scales ===============
% For M0/A, use the SMALLEST positive denominator among the three test
% years at each grid cell; for B, the same fixed 2004-2020 historical IQR
% serves all years. This explicitly targets small denominators, without
% using model errors to select cells.
requiredStage1 = { ...
    'comparisonMask','denom_M0_annual','denom_A_annual','denom_B_fixed'};
checkMATVariables(stage1File,requiredStage1);
D = load(stage1File,requiredStage1{:});
assert(isequal(logical(D.comparisonMask),comparisonMask), ...
    'Stage 1 and Stage 2 comparison masks are not identical.');
assert(isequal(size(D.denom_M0_annual),[numel(lat),numel(lon),3]), ...
    'Unexpected Stage 1 M0 annual denominator array dimensions.');
assert(isequal(size(D.denom_A_annual),[numel(lat),numel(lon),3]), ...
    'Unexpected Stage 1 A annual denominator array dimensions.');
assert(isequal(size(D.denom_B_fixed),[numel(lat),numel(lon)]), ...
    'Unexpected Stage 1 fixed historical-IQR array dimensions.');

smallDenomMaps = { ...
    min(double(D.denom_M0_annual),[],3); ...
    min(double(D.denom_A_annual),[],3); ...
    double(D.denom_B_fixed)};
clear D;
for k=1:nScheme
    dd = smallDenomMaps{k};
    if ~all(isfinite(dd(comparisonMask)) & dd(comparisonMask)>0)
        error('Nonfinite or nonpositive denominators for %s.',schemeNames{k});
    end
end

%% =============== 4. Reconstruct original components ===================
% The entire original six-class decomposition is reproduced separately
% for M0/A/B. No per-model predictions are loaded or altered here.
A = cell(nScheme,1);
Table1_all = cell(nScheme,1);
classStatistics_all = cell(nScheme,1);

for k=1:nScheme
    A{k} = reconstructComponents(R_input{k},AR_input{k}, ...
        comparisonMask,classNames,schemeNames{k});
    Table1_all{k} = A{k}.Table1;
    classStatistics_all{k} = A{k}.ClassStatistics;
    if ~isequal(A{k}.PixelCount,originalCategoryCounts(:,k))
        error('Stage 3 class counts do not match Stage 2 for %s.', ...
            schemeNames{k});
    end
    if A{k}.nAssigned~=nValid || A{k}.nUnassigned~=0 || ...
            A{k}.nMultiAssigned~=0
        error('Invalid category assignments under %s.',schemeNames{k});
    end
    if A{k}.maxAbsReconstructionError>1e-8
        error('Gain decomposition fails telescoping identity for %s.', ...
            schemeNames{k});
    end
    fprintf('%-23s: assigned=%5d; max reconstruction error=%.3g\n', ...
        schemeNames{k},A{k}.nAssigned,A{k}.maxAbsReconstructionError);
end

A_M0 = A{1}; A_A = A{2}; A_B = A{3};
Table1_M0 = Table1_all{1};
Table1_A = Table1_all{2};
Table1_B = Table1_all{3};
ClassStatistics_M0 = classStatistics_all{1};
ClassStatistics_A = classStatistics_all{2};
ClassStatistics_B = classStatistics_all{3};

%% =============== 5. Audit M0 against manuscript Table 1 ==============
% This is a strict reproduction audit, not the sensitivity comparison.
% Class 1 (No gain) has zero total reduction and undefined shares.
auditClass = classNames;
auditPixelCount = A_M0.PixelCount;
maxShareDiffPP = nan(nClass,1);
maxMeanDiff = nan(nClass,1);
auditStatus = cell(nClass,1);

for c=1:nClass
    observedShare = [ ...
        Table1_M0.LinearContributionPct(c), ...
        Table1_M0.NonlinearContributionPct(c), ...
        Table1_M0.HydroContributionPct(c), ...
        Table1_M0.MixedContributionPct(c)];
    observedMean = [ ...
        Table1_M0.MeanLinearReduction(c), ...
        Table1_M0.MeanNonlinearReduction(c), ...
        Table1_M0.MeanHydroReduction(c), ...
        Table1_M0.MeanMixedReduction(c)];
    expectedShare = expectedM0Shares(c,:);
    expectedMean = expectedM0MeanComponents(c,:);
    if c==1
        shareOK = all(isnan(observedShare));
        maxShareDiffPP(c) = 0;
    else
        maxShareDiffPP(c) = max(abs(observedShare-expectedShare));
        shareOK = maxShareDiffPP(c)<=auditShareTolPP;
    end
    maxMeanDiff(c) = max(abs(observedMean-expectedMean));
    countOK = auditPixelCount(c)==expectedM0Counts(c);
    meanOK = maxMeanDiff(c)<=auditMeanTol;
    if countOK && shareOK && meanOK
        auditStatus{c} = 'PASS';
    else
        auditStatus{c} = 'FAIL';
    end
end

M0AuditTable = table(auditClass,auditPixelCount,expectedM0Counts, ...
    maxShareDiffPP,maxMeanDiff,auditStatus, ...
    'VariableNames',{'Category','ComputedCount','ManuscriptCount', ...
    'MaxShareDifference_pp','MaxMeanReductionDifference','Status'});
if ~all(strcmp(auditStatus,'PASS'))
    disp(M0AuditTable);
    error(['M0 does not reproduce the original manuscript Table 1. ', ...
        'Please resolve the audit before interpreting A/B results.']);
end
fprintf('M0 MANUSCRIPT TABLE 1 AUDIT: 6/6 PASS.\n\n');

%% =========== 6. Main hydroclimatic-category comparison ===============
% History contribution = Linear + Nonlinear WITHIN category 4.
hydroMain = table(schemeNames,zeros(nScheme,1), ...
    nan(nScheme,1),nan(nScheme,1),nan(nScheme,1), ...
    nan(nScheme,1),nan(nScheme,1),nan(nScheme,1), ...
    nan(nScheme,1), ...
    'VariableNames',{'Scheme','HydroClassPixels', ...
    'LinearPct','NonlinearPct','HistoryPct', ...
    'HydroclimaticPct','MixedPct','MeanTotalReduction', ...
    'HydroMinusM0_pp'});

for k=1:nScheme
    T=Table1_all{k};
    c=hydroClassIdx;
    hydroMain.HydroClassPixels(k)=A{k}.PixelCount(c);
    hydroMain.LinearPct(k)=T.LinearContributionPct(c);
    hydroMain.NonlinearPct(k)=T.NonlinearContributionPct(c);
    hydroMain.HistoryPct(k)=T.LinearContributionPct(c)+T.NonlinearContributionPct(c);
    hydroMain.HydroclimaticPct(k)=T.HydroContributionPct(c);
    hydroMain.MixedPct(k)=T.MixedContributionPct(c);
    hydroMain.MeanTotalReduction(k)=A{k}.MeanTotalReduction(c);
end
hydroMain.HydroMinusM0_pp=hydroMain.HydroclimaticPct-hydroMain.HydroclimaticPct(1);

%% =================== 7. Component share tables =======================
% One row per scheme/category; each ratio is the ratio of component sums.
longScheme = cell(nScheme*nClass,1);
longCategory = cell(nScheme*nClass,1);
longCount = zeros(nScheme*nClass,1);
longLinear = nan(nScheme*nClass,1);
longNonlinear = nan(nScheme*nClass,1);
longHydro = nan(nScheme*nClass,1);
longMixed = nan(nScheme*nClass,1);
longTotal = nan(nScheme*nClass,1);
longMeanTotal = nan(nScheme*nClass,1);
longMeanRelTotal = nan(nScheme*nClass,1);
row = 0;
for k=1:nScheme
    T=Table1_all{k};
    CS=classStatistics_all{k};
    for c=1:nClass
        row=row+1;
        longScheme{row}=schemeNames{k};
        longCategory{row}=classNames{c};
        longCount(row)=A{k}.PixelCount(c);
        longLinear(row)=T.LinearContributionPct(c);
        longNonlinear(row)=T.NonlinearContributionPct(c);
        longHydro(row)=T.HydroContributionPct(c);
        longMixed(row)=T.MixedContributionPct(c);
        longTotal(row)=CS.TotalReductionSum(c);
        longMeanTotal(row)=A{k}.MeanTotalReduction(c);
        longMeanRelTotal(row)=A{k}.MeanTotalReductionPct(c);
    end
end
componentShareComparison = table(longScheme,longCategory,longCount, ...
    longLinear,longNonlinear,longHydro,longMixed,longTotal, ...
    'VariableNames',{'Scheme','Category','Pixels','LinearPct', ...
    'NonlinearPct','HydroclimaticPct','MixedPct','TotalReductionSum'});
categoryStatisticsComparison = table(longScheme,longCategory,longCount, ...
    longMeanTotal,longMeanRelTotal,longTotal, ...
    'VariableNames',{'Scheme','Category','Pixels', ...
    'MeanTotalReduction','MeanRelativeTotalReductionPct', ...
    'TotalReductionSum'});

%% =============== 8. Extreme-cell influence diagnostics ================
% Diagnostics DO NOT alter main Table 1 or saved main component maps.
% Main shares pool ALL cells; exclusions below are transparent, secondary
% tests applied separately to each scheme and each positive-gain class.
%   Full                  : original full class
%   ExcludeMaximumGain    : remove the single cell with largest total gain
%   ExcludeTop1PctGain    : remove ceil(1% * class size) largest-gain cells
%   ExcludeSmallest1PctDenom: remove ceil(1% * class size) cells with
%                            smallest denominators, independent of error
% Small-denominator ranking:
%   M0: smallest 2022/2023/2024 test-year Range at that grid cell
%   A : smallest 2022/2023/2024 test-year IQR at that grid cell
%   B : fixed 2004-2020 historical IQR

numRows = nScheme*(nClass-1)*4;
diagScheme = cell(numRows,1);
diagCategory = cell(numRows,1);
diagMethod = cell(numRows,1);
diagOriginalCount = zeros(numRows,1);
diagRemoved = zeros(numRows,1);
diagRemaining = zeros(numRows,1);
diagTotalGain = nan(numRows,1);
diagLinearShare = nan(numRows,1);
diagNonlinearShare = nan(numRows,1);
diagHydroShare = nan(numRows,1);
diagHistoryShare = nan(numRows,1);
diagMixedShare = nan(numRows,1);
diagHydroChangeFromFull = nan(numRows,1);

concRows = nScheme*(nClass-1);
concScheme = cell(concRows,1);
concClass = cell(concRows,1);
concPixels = zeros(concRows,1);
concTop1Count = zeros(concRows,1);
concTop1GainShare = nan(concRows,1);
concMaxCellGainShare = nan(concRows,1);
concSmallDenom1GainShare = nan(concRows,1);
concMinDenom = nan(concRows,1);

methods = {'Full';'ExcludeMaximumGain';'ExcludeTop1PctGain'; ...
           'ExcludeSmallest1PctDenom'};
r=0;
q=0;
for k=1:nScheme
    for c=2:nClass
        thisClassMask = A{k}.classMasks{c};
        totalGain = A{k}.G_total(thisClassMask);
        linearGain = A{k}.G_linear(thisClassMask);
        nonlinearGain = A{k}.G_nonlinear(thisClassMask);
        hydroGain = A{k}.G_hydro(thisClassMask);
        mixedGain = A{k}.G_mixed(thisClassMask);
        denomMap = smallDenomMaps{k};
        denominator = denomMap(thisClassMask);
        n=numel(totalGain);
        if n==0
            error('No cells in positive-gain class %s under %s.', ...
                classNames{c},schemeNames{k});
        end
        [~,descendingGainOrder]=sort(totalGain,'descend');
        [~,ascendingDenomOrder]=sort(denominator,'ascend');
        nTop=max(1,ceil(extremeFraction*n));
        keepMasks=cell(4,1);
        keepMasks{1}=true(n,1);
        keepMasks{2}=true(n,1);
        keepMasks{2}(descendingGainOrder(1))=false;
        keepMasks{3}=true(n,1);
        keepMasks{3}(descendingGainOrder(1:nTop))=false;
        keepMasks{4}=true(n,1);
        keepMasks{4}(ascendingDenomOrder(1:nTop))=false;

        fullShare = componentShares([ ...
            sum(linearGain),sum(nonlinearGain), ...
            sum(hydroGain),sum(mixedGain)],sum(totalGain));
        q=q+1;
        concScheme{q}=schemeNames{k};
        concClass{q}=classNames{c};
        concPixels(q)=n;
        concTop1Count(q)=nTop;
        concTop1GainShare(q)=sum(totalGain(descendingGainOrder(1:nTop))) ...
            /sum(totalGain)*100;
        concMaxCellGainShare(q)=totalGain(descendingGainOrder(1)) ...
            /sum(totalGain)*100;
        concSmallDenom1GainShare(q)=sum(totalGain(ascendingDenomOrder(1:nTop))) ...
            /sum(totalGain)*100;
        concMinDenom(q)=min(denominator);

        for j=1:numel(methods)
            keep=keepMasks{j};
            componentSums = [sum(linearGain(keep)), ...
                sum(nonlinearGain(keep)),sum(hydroGain(keep)), ...
                sum(mixedGain(keep))];
            tot = sum(totalGain(keep));
            shares=componentShares(componentSums,tot);
            r=r+1;
            diagScheme{r}=schemeNames{k};
            diagCategory{r}=classNames{c};
            diagMethod{r}=methods{j};
            diagOriginalCount(r)=n;
            diagRemoved(r)=nnz(~keep);
            diagRemaining(r)=nnz(keep);
            diagTotalGain(r)=tot;
            diagLinearShare(r)=shares(1);
            diagNonlinearShare(r)=shares(2);
            diagHydroShare(r)=shares(3);
            diagHistoryShare(r)=shares(1)+shares(2);
            diagMixedShare(r)=shares(4);
            diagHydroChangeFromFull(r)=shares(3)-fullShare(3);
        end
    end
end

extremeDiagnostics = table(diagScheme,diagCategory,diagMethod, ...
    diagOriginalCount,diagRemoved,diagRemaining,diagTotalGain, ...
    diagLinearShare,diagNonlinearShare,diagHistoryShare,diagHydroShare, ...
    diagMixedShare,diagHydroChangeFromFull, ...
    'VariableNames',{'Scheme','Category','Method','OriginalPixels', ...
    'ExcludedPixels','RetainedPixels','RetainedTotalReductionSum', ...
    'LinearPct','NonlinearPct','HistoryPct', ...
    'HydroclimaticPct','MixedPct','HydroMinusFull_pp'});

extremeConcentration = table(concScheme,concClass,concPixels, ...
    concTop1Count,concTop1GainShare,concMaxCellGainShare, ...
    concSmallDenom1GainShare,concMinDenom, ...
    'VariableNames',{'Scheme','Category','ClassPixels','Top1PctPixels', ...
    'Top1PctGainShareOfTotalPct','LargestGainCellSharePct', ...
    'Smallest1PctDenomCellsGainSharePct','MinimumDenominator_mm'});

% Concise hydroclimatic-only diagnostic: directly reveals whether the
% pooled 71%/29% is dominated by a small fraction of high-gain cells.
hydroExtremeRows = strcmp(extremeDiagnostics.Category, ...
    classNames{hydroClassIdx});
hydroExtremeDiagnostics = extremeDiagnostics(hydroExtremeRows,:);
hydroConcentrationRows = strcmp(extremeConcentration.Category, ...
    classNames{hydroClassIdx});
hydroExtremeConcentration = extremeConcentration(hydroConcentrationRows,:);

%% ======================== 9. Display results ==========================
fprintf('\n=============== M0 TABLE 1 (REPRODUCTION) ===============\n');
disp(Table1_M0);
fprintf('\n=============== A TABLE 1 (TEST IQR) ====================\n');
disp(Table1_A);
fprintf('\n=============== B TABLE 1 (HISTORY IQR) =================\n');
disp(Table1_B);
fprintf('\n=============== HYDRO-CLASS CONTRIBUTIONS =================\n');
disp(hydroMain);
fprintf('\n=============== HYDRO-CLASS CONCENTRATION ================\n');
disp(hydroExtremeConcentration);
fprintf('\n=============== HYDRO-CLASS EXTREME DIAGNOSTICS =========\n');
disp(hydroExtremeDiagnostics);
fprintf('\nM0 TABLE 1 AUDIT: %d/%d PASS.\n', ...
    nnz(strcmp(auditStatus,'PASS')),nClass);

%% ========================== 10. Save ================================
matFile=fullfile(outDir,'s5_NRMSE_sensitivity3_results.mat');
save(matFile, ...
    'A_M0','A_A','A_B', ...
    'Table1_M0','Table1_A','Table1_B', ...
    'ClassStatistics_M0','ClassStatistics_A','ClassStatistics_B', ...
    'hydroMain','componentShareComparison','categoryStatisticsComparison', ...
    'extremeDiagnostics','extremeConcentration', ...
    'hydroExtremeDiagnostics','hydroExtremeConcentration', ...
    'M0AuditTable','comparisonMask','modelNames', ...
    'schemeNames','classNames','componentNames', ...
    'testYears','lat','lon','stage1File','stage2File','extremeFraction', ...
    '-v7.3');

writetable(Table1_M0,fullfile(outDir,'s5_stage3_Table1_M0.csv'));
writetable(Table1_A,fullfile(outDir,'s5_stage3_Table1_A.csv'));
writetable(Table1_B,fullfile(outDir,'s5_stage3_Table1_B.csv'));
writetable(hydroMain,fullfile(outDir,'s5_stage3_Hydro_main_comparison.csv'));
writetable(componentShareComparison, ...
    fullfile(outDir,'s5_stage3_component_shares_all_categories.csv'));
writetable(categoryStatisticsComparison, ...
    fullfile(outDir,'s5_stage3_category_statistics.csv'));
writetable(extremeDiagnostics, ...
    fullfile(outDir,'s5_stage3_extreme_cell_diagnostics.csv'));
writetable(extremeConcentration, ...
    fullfile(outDir,'s5_stage3_extreme_cell_concentration.csv'));
writetable(M0AuditTable,fullfile(outDir,'s5_stage3_M0_audit.csv'));

excelFile=fullfile(outDir,'s5_NRMSE_sensitivity3_tables.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(Table1_M0,excelFile,'Sheet','Table1_M0');
writetable(Table1_A,excelFile,'Sheet','Table1_A');
writetable(Table1_B,excelFile,'Sheet','Table1_B');
writetable(hydroMain,excelFile,'Sheet','Hydro_Main');
writetable(componentShareComparison,excelFile,'Sheet','AllClass_Shares');
writetable(categoryStatisticsComparison,excelFile,'Sheet','Class_Statistics');
writetable(hydroExtremeDiagnostics,excelFile,'Sheet','Hydro_Extreme');
writetable(hydroExtremeConcentration,excelFile,'Sheet','Hydro_Concentration');
writetable(extremeDiagnostics,excelFile,'Sheet','All_Extreme');
writetable(extremeConcentration,excelFile,'Sheet','All_Concentration');
writetable(M0AuditTable,excelFile,'Sheet','M0_Audit');

if makeFigure
    try
        % Main FULL-DOMAIN component shares only (no diagnostic exclusions).
        figureHandle=figure('Color','w','Position',[100,100,960,540], ...
            'Visible','off');
        seriesValues=[hydroMain.LinearPct,hydroMain.NonlinearPct, ...
            hydroMain.HydroclimaticPct];
        bars=bar(seriesValues,'stacked','BarWidth',0.55);
        bars(1).FaceColor=[0.75,0.27,0.21];
        bars(2).FaceColor=[0.20,0.43,0.73];
        bars(3).FaceColor=[0.18,0.55,0.34];
        ax=gca;
        ax.XTick=1:nScheme;
        ax.XTickLabel=shortSchemeNames;
        ax.FontName='Arial';
        ax.FontSize=12;
        ylim([0,105]);
        ylabel('Share of total NRMSE reduction (%)');
        xlabel('Normalization scheme');
        legend({'Linear','Nonlinear','Hydroclimatic'}, ...
            'Location','eastoutside','Box','off');
        title(['Hydroclimatic-gain category: ', ...
            'winner-conditioned component contributions'], ...
            'FontWeight','normal');
        grid on;
        ax.GridAlpha=0.12;
        exportgraphics(figureHandle,fullfile(outDir, ...
            'Fig1_Hydro_Contribution_Sensitivity.png'),'Resolution',300);
        close(figureHandle);
    catch ME
        warning('Optional hydro-component figure could not be created: %s', ...
            ME.message);
        if exist('figureHandle','var') && isgraphics(figureHandle)
            close(figureHandle);
        end
    end
end

fprintf('\nSaved MAT: %s\n',matFile);
fprintf('Saved Excel: %s\n',excelFile);
fprintf('Elapsed: %.1f s\n',toc(startTic));
fprintf('STAGE 3 COMPLETED: Table 1 sensitivity and influence diagnostics.\n');
fprintf(['MAIN TABLES ALWAYS INCLUDE EVERY VALID CELL. ', ...
    'EXCLUSIONS ARE DIAGNOSTIC ONLY.\n']);
fprintf('===============================================================\n');
diary off;

%% ======================== Local functions ============================
function checkMATVariables(filename,variableNames)
    names=who('-file',filename);
    missing=setdiff(variableNames,names);
    if ~isempty(missing)
        error('MAT file %s is missing: %s', ...
            filename,strjoin(missing,', '));
    end
end

function shares = componentShares(componentSums,totalGainSum)
    % Ratio of component SUMS to total-gain SUM, not the mean of ratios.
    if ~isfinite(totalGainSum) || abs(totalGainSum)<=1e-12
        shares=nan(1,4);
    else
        shares=double(componentSums)./totalGainSum.*100;
    end
end

function A = reconstructComponents(R,AR,analysisMask,classNames,schemeName)
    % Exact reproduction of original s2_gain_attribution3.m (Original run).
    % Class order: No gain, Linear, Nonlinear, Hydro, Mixed-L, Mixed-N.
    assert(isequal(logical(R.analysisMask),analysisMask), ...
        'Analysis-mask mismatch for %s.',schemeName);
    assert(isequal(size(AR),size(analysisMask)), ...
        'AR map does not match Stage 2 mask for %s.',schemeName);

    L=double(R.linearFusion);
    N=double(R.nonlinearFusion);
    X=double(R.nonlinearXFusion);
    B=double(R.allFusion);
    assert(isequal(size(L),size(analysisMask)) && ...
        isequal(size(N),size(analysisMask)) && ...
        isequal(size(X),size(analysisMask)) && ...
        isequal(size(B),size(analysisMask)), ...
        'Fusion arrays have incorrect dimensions for %s.',schemeName);

    if ~all(isfinite(AR(analysisMask)) & isfinite(L(analysisMask)) & ...
            isfinite(N(analysisMask)) & isfinite(X(analysisMask)) & ...
            isfinite(B(analysisMask)))
        error('Nonfinite model performance within shared mask for %s.',schemeName);
    end

    classMasks={ ...
        analysisMask & logical(R.maskNoGain); ...
        analysisMask & logical(R.maskLinearGain); ...
        analysisMask & logical(R.maskNonlinearGain); ...
        analysisMask & logical(R.maskHydroClean); ...
        analysisMask & logical(R.maskMixedLinearSecond); ...
        analysisMask & logical(R.maskMixedNonlinearSecond)};

    nClass=numel(classNames);
    classSum=zeros(size(analysisMask));
    for c=1:nClass
        classSum=classSum+double(classMasks{c});
        if nnz(classMasks{c} & R.gainClass6Map~=c)>0
            error('Stage 2 class map disagrees with mask %d under %s.', ...
                c,schemeName);
        end
    end
    nAssigned=nnz(analysisMask & classSum==1);
    nUnassigned=nnz(analysisMask & classSum==0);
    nMultiAssigned=nnz(analysisMask & classSum>1);
    if nUnassigned>0 || nMultiAssigned>0
        error('Overlapping or missing class assignments for %s.',schemeName);
    end

    G_total=AR-B;
    G_total(~analysisMask)=NaN;

    G_linear=nan(size(AR));
    G_nonlinear=nan(size(AR));
    G_hydro=nan(size(AR));
    G_mixed=nan(size(AR));
    G_linear(analysisMask)=0;
    G_nonlinear(analysisMask)=0;
    G_hydro(analysisMask)=0;
    G_mixed(analysisMask)=0;

    % Class 1: No gain (AR wins), all components remain zero.
    % Class 2: Linear gain.
    q=classMasks{2};
    G_linear(q)=AR(q)-B(q);

    % Class 3: Nonlinear gain, telescoping AR -> best Linear -> Nonlinear.
    q=classMasks{3};
    G_linear(q)=AR(q)-L(q);
    G_nonlinear(q)=L(q)-N(q);

    % Class 4: Hydroclimatic gain, AR -> L -> N -> X (paired winner).
    q=classMasks{4};
    G_linear(q)=AR(q)-L(q);
    G_nonlinear(q)=L(q)-N(q);
    G_hydro(q)=N(q)-X(q);

    % Class 5: Mixed-L, AR -> L -> X.
    q=classMasks{5};
    G_linear(q)=AR(q)-L(q);
    G_mixed(q)=L(q)-X(q);

    % Class 6: Mixed-N, AR -> L -> N -> X (unpaired winner).
    q=classMasks{6};
    G_linear(q)=AR(q)-L(q);
    G_nonlinear(q)=L(q)-N(q);
    G_mixed(q)=N(q)-X(q);

    % Original s2_gain_attribution3.m removes only tiny roundoff noise.
    tolerance=1e-10;
    G_total(abs(G_total)<tolerance)=0;
    G_linear(abs(G_linear)<tolerance)=0;
    G_nonlinear(abs(G_nonlinear)<tolerance)=0;
    G_hydro(abs(G_hydro)<tolerance)=0;
    G_mixed(abs(G_mixed)<tolerance)=0;

    recon=G_linear+G_nonlinear+G_hydro+G_mixed;
    diff=recon-G_total;
    maxAbsReconstructionError=max(abs(diff(analysisMask)));
    if any(G_total(analysisMask)<-1e-8)
        error('Negative total gain beyond rounding error under %s.',schemeName);
    end

    PixelCount=zeros(nClass,1);
    PixelPct=zeros(nClass,1);
    MeanTotalReduction=nan(nClass,1);
    MeanTotalReductionPct=nan(nClass,1);
    MeanLinearReduction=nan(nClass,1);
    MeanNonlinearReduction=nan(nClass,1);
    MeanHydroReduction=nan(nClass,1);
    MeanMixedReduction=nan(nClass,1);
    LinearContributionPct=nan(nClass,1);
    NonlinearContributionPct=nan(nClass,1);
    HydroContributionPct=nan(nClass,1);
    MixedContributionPct=nan(nClass,1);
    TotalReductionSum=zeros(nClass,1);

    nValid=nnz(analysisMask);
    for c=1:nClass
        q=classMasks{c};
        PixelCount(c)=nnz(q);
        PixelPct(c)=PixelCount(c)/nValid*100;
        if PixelCount(c)==0
            continue;
        end
        tot=G_total(q);
        l=G_linear(q);
        n=G_nonlinear(q);
        h=G_hydro(q);
        m=G_mixed(q);
        MeanTotalReduction(c)=mean(tot,'omitnan');
        validRelativeMask=q & isfinite(AR) & AR>0;
        if nnz(validRelativeMask)>0
            MeanTotalReductionPct(c)=mean( ...
                G_total(validRelativeMask)./AR(validRelativeMask).*100, ...
                'omitnan');
        end
        MeanLinearReduction(c)=mean(l,'omitnan');
        MeanNonlinearReduction(c)=mean(n,'omitnan');
        MeanHydroReduction(c)=mean(h,'omitnan');
        MeanMixedReduction(c)=mean(m,'omitnan');
        TotalReductionSum(c)=sum(tot,'omitnan');
        shares=componentShares([sum(l,'omitnan'),sum(n,'omitnan'), ...
            sum(h,'omitnan'),sum(m,'omitnan')],TotalReductionSum(c));
        LinearContributionPct(c)=shares(1);
        NonlinearContributionPct(c)=shares(2);
        HydroContributionPct(c)=shares(3);
        MixedContributionPct(c)=shares(4);
    end

    Table1=table(classNames,LinearContributionPct,NonlinearContributionPct, ...
        HydroContributionPct,MixedContributionPct, ...
        MeanLinearReduction,MeanNonlinearReduction, ...
        MeanHydroReduction,MeanMixedReduction, ...
        'VariableNames',{'AttributionClass','LinearContributionPct', ...
        'NonlinearContributionPct','HydroContributionPct', ...
        'MixedContributionPct','MeanLinearReduction', ...
        'MeanNonlinearReduction','MeanHydroReduction', ...
        'MeanMixedReduction'});

    ClassStatistics=table(classNames,PixelCount,PixelPct, ...
        MeanTotalReduction,MeanTotalReductionPct,TotalReductionSum, ...
        'VariableNames',{'AttributionClass','Pixels','PixelPct', ...
        'MeanTotalReduction','MeanTotalReductionPct','TotalReductionSum'});

    A.schemeName=schemeName;
    A.analysisMask=analysisMask;
    A.nAssigned=nAssigned;
    A.nUnassigned=nUnassigned;
    A.nMultiAssigned=nMultiAssigned;
    A.nValid=nValid;
    A.maxAbsReconstructionError=maxAbsReconstructionError;
    A.classMasks=classMasks;
    A.PixelCount=PixelCount;
    A.MeanTotalReduction=MeanTotalReduction;
    A.MeanTotalReductionPct=MeanTotalReductionPct;
    A.AR_NRMSE=AR;
    A.G_total=G_total;
    A.G_linear=G_linear;
    A.G_nonlinear=G_nonlinear;
    A.G_hydro=G_hydro;
    A.G_mixed=G_mixed;
    A.Table1=Table1;
    A.ClassStatistics=ClassStatistics;
end
