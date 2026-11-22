%% s5_NRMSE_sensitivity2.m
% =========================================================================
% Stage 2: Winner and gain-CATEGORY sensitivity to NRMSE normalization.
%
% Input:
%   s5_NRMSE_sensitivity1_results/s5_NRMSE_sensitivity1_results.mat
%   (Original-run nine-model results; NO forecasting/training is repeated.)
%
% Three NRMSE definitions carried over UNCHANGED from Stage 1:
%   M0: annual test RMSE / that test year's observed Range
%   A : annual test RMSE / that test year's observed IQR
%   B : annual test RMSE / FIXED 2004-2020 observed IQR
%
% IMPORTANT: All nine models are compared using each model's THREE-YEAR
% mean of independently calculated annual NRMSE, NOT a pooled 36-month RMSE.
% The same 17,284-pixel comparisonMask is used under M0, A, and B.
%
% Categories reproduce original scripts s2_gain_attribution2.m (five
% combined groups) and s2_gain_attribution3.m (six detailed categories):
%   1 No gain       : Linear best, and within Linear AR is best.
%   2 Linear gain   : Linear best, and ARIMA or SARIMA is best.
%   3 Nonlinear gain: Nonlinear group is best (SVM, RF, LSTM).
%   4 Hydro gain    : Hydro-X group best, Nonlinear group second-best,
%                    AND winning Hydro-X and Nonlinear algorithms paired.
%   5 Mixed-L       : Hydro-X group best, Linear group second-best.
%   6 Mixed-N       : Hydro-X best, Nonlinear second-best, NOT paired.
%
% A direct pairing means SVM<->SVMX, RF<->RFX, or LSTM<->LSTMX.
% The winning model and group tie rule is the original MATLAB min/sort
% rule (first member of the defined order for an exact tie).
%
% OUTPUTS (new directory; all existing results are left untouched):
%   s5_NRMSE_sensitivity2_results.mat   : maps and component input structs
%   s5_NRMSE_sensitivity2_tables.xlsx  : composition, agreement, transition
%   s5_stage2_winner_composition.csv
%   s5_stage2_category_composition.csv
%   s5_stage2_group_composition.csv
%   s5_stage2_agreement.csv
%   s5_stage2_diagnostics.csv
%   s5_stage2_M0_audit.csv
%   s5_stage2_transitions_*.csv        : counts and row percentages
%   s5_NRMSE_sensitivity2_log.txt
%
% NOT in this stage: attribution-component shares (e.g. 71%/29%), which
% will be recalculated in Stage 3 using saved R_M0, R_A, R_B structs.
% =========================================================================

clc;
clearvars;

%% ================== 1. User settings and output paths ==================
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
baseDir = fullfile(projectRoot, 's3_darft', 's1_analysis');
stage1File = fullfile(baseDir, 's5_NRMSE_sensitivity1_results', ...
    's5_NRMSE_sensitivity1_results.mat');

outDir = fullfile(baseDir, 's5_NRMSE_sensitivity2_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% M0 diagnostic only: try reading the original results of s2 scripts.
% These historical MAT files are NOT needed for calculating A/B.
oldFiveClassFile = fullfile(baseDir, 's2_gain_attribution2_results', ...
    's2_gain_attribution2_results.mat');
oldSixClassFile = fullfile(baseDir, 's2_gain_attribution3_results', ...
    's2_gain_attribution3_results.mat');

expectedValidCount = 17284;    % Warning only; original comparisonMask rules.
expectedM0ClassPct = [8.1;25.1;26.1;12.1;17.4;11.2];
expectedM0ClassTolerancePP = 0.055;  % One-decimal manuscript rounding.

modelNamesExpected = { ...
    'AR';'ARIMA';'SARIMA';'SVM';'RF';'LSTM';'SVMX';'RFX';'LSTMX'};
schemeNamesExpected = { ...
    'M0_test_range';'A_test_IQR';'B_fixed_history_IQR'};
categoryNames = { ...
    'No gain';'Linear gain';'Nonlinear gain';'Hydroclimatic gain'; ...
    'Mixed-L';'Mixed-N'};
modelGroupNames = {'Linear';'Nonlinear';'Hydroclimatic-X'};
pairNames = {'M0_vs_A';'M0_vs_B';'A_vs_B'};
pairSchemeIdx = [1 2;1 3;2 3];

nModel = numel(modelNamesExpected);
nClass = numel(categoryNames);
nScheme = numel(schemeNamesExpected);
nPair = size(pairSchemeIdx,1);

logFile = fullfile(outDir,'s5_NRMSE_sensitivity2_log.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
startTic = tic;

fprintf('=============================================================\n');
fprintf('S5 STAGE 2: WINNER AND GAIN-CATEGORY SENSITIVITY\n');
fprintf('=============================================================\n');
fprintf('M0 = test-year Range; A = test-year IQR; B = fixed 2004-2020 IQR\n');
fprintf('Stage 1 input: %s\n\n',stage1File);

%% ======================== 2. Load Stage 1 ============================
if ~isfile(stage1File)
    error('Stage 1 MAT file is missing:\n%s',stage1File);
end

needed = { ...
    'modelNames','schemeNames','testYears','lat','lon', ...
    'comparisonMask','M0ReferenceMask', ...
    'NRMSE_mean_M0','NRMSE_mean_A','NRMSE_mean_B', ...
    'M0AuditTable'};
existing = who('-file',stage1File);
missing = setdiff(needed,existing);
if ~isempty(missing)
    error('Stage 1 MAT is missing variable(s): %s',strjoin(missing,', '));
end

S = load(stage1File,needed{:});
modelNames = S.modelNames(:);
schemeNames = S.schemeNames(:);
testYears = double(S.testYears(:)');
lat = double(S.lat(:));
lon = double(S.lon(:));
comparisonMask = logical(S.comparisonMask);
M0ReferenceMask = logical(S.M0ReferenceMask);
M0AuditTable = S.M0AuditTable;

assert(isequal(modelNames,modelNamesExpected), ...
    'Model order is not AR, ARIMA, SARIMA, SVM, RF, LSTM, SVMX, RFX, LSTMX.');
assert(isequal(schemeNames,schemeNamesExpected), ...
    'Scheme names/order do not match the Stage 1 M0/A/B definitions.');
assert(isequal(testYears,[2022 2023 2024]), ...
    'Unexpected test-year sequence.');

% Never silently continue when original evaluation could not be reproduced.
assert(istable(M0AuditTable) && ...
       ismember('AuditStatus',M0AuditTable.Properties.VariableNames), ...
    'Stage 1 M0AuditTable/AuditStatus is unavailable.');
if ~all(strcmp(M0AuditTable.AuditStatus, 'PASS'))
    error('Stage 1 M0 reproduction failed. Resolve its audit before Stage 2.');
end

[nLat,nLon] = size(comparisonMask);
assert(nLat == numel(lat) && nLon == numel(lon), ...
    'Comparison mask spatial size does not match lat/lon.');
nPixel = nnz(comparisonMask);
if nPixel ~= expectedValidCount
    warning('Stage 2 common mask has %d cells; expected %d.', ...
        nPixel,expectedValidCount);
end
if any(comparisonMask(:) & ~M0ReferenceMask(:))
    error('comparisonMask contains cells outside the Stage 1 M0 domain.');
end

scoreMaps = cell(nScheme,1);
scoreMaps{1} = double(S.NRMSE_mean_M0);
scoreMaps{2} = double(S.NRMSE_mean_A);
scoreMaps{3} = double(S.NRMSE_mean_B);
clear S;
for k = 1:nScheme
    if ~isequal(size(scoreMaps{k}),[nLat,nLon,nModel])
        error('Mean NRMSE array for %s has unexpected dimensions.', ...
            schemeNames{k});
    end
    for m = 1:nModel
        x = scoreMaps{k}(:,:,m);
        if ~all(isfinite(x(comparisonMask)))
            error('Nonfinite NRMSE in shared mask: %s, %s.', ...
                schemeNames{k},modelNames{m});
        end
    end
end

fprintf('Common valid grid cells : %d\n',nPixel);
fprintf('Test years              : %s\n',mat2str(testYears));
fprintf('Stage 1 M0 audit        : 9/9 PASS\n\n');

%% ================= 3. Winner and classification maps ==================
% R_* fields deliberately match the names used by the original s2 scripts,
% so Stage 3 can reconstruct components without changing decomposition.
R = cell(nScheme,1);
for k = 1:nScheme
    R{k} = classifyThreeYearWinners(scoreMaps{k},comparisonMask);
    if R{k}.nAssigned ~= nPixel || R{k}.nUnassigned ~= 0 || ...
            R{k}.nMultiAssigned ~= 0 || R{k}.nWinnerMismatch ~= 0
        error(['Category/winner reconstruction invalid in %s: ', ...
            'assigned=%d unassigned=%d multi=%d winner-mismatch=%d.'], ...
            schemeNames{k},R{k}.nAssigned,R{k}.nUnassigned, ...
            R{k}.nMultiAssigned,R{k}.nWinnerMismatch);
    end
    fprintf('%-25s | classified %5d | unassigned %d | model mismatch %d\n', ...
        schemeNames{k},R{k}.nAssigned,R{k}.nUnassigned, ...
        R{k}.nWinnerMismatch);
end
fprintf('\n');

% The original scripts store a 5-category map (Mixed-L/N combined). Stage
% 2 stores BOTH the 5-category map and the six detailed categories.
R_M0 = R{1};
R_A = R{2};
R_B = R{3};

%% ================ 4. Per-scheme composition counts ====================
winnerCounts = zeros(nModel,nScheme);
winnerPct = zeros(nModel,nScheme);
categoryCounts = zeros(nClass,nScheme);
categoryPct = zeros(nClass,nScheme);
groupCounts = zeros(numel(modelGroupNames),nScheme);
groupPct = zeros(numel(modelGroupNames),nScheme);

for k = 1:nScheme
    for m = 1:nModel
        winnerCounts(m,k) = nnz(comparisonMask & R{k}.bestAll9Idx==m);
    end
    for c = 1:nClass
        categoryCounts(c,k) = nnz(comparisonMask & R{k}.gainClass6Map==c);
    end
    for g = 1:numel(modelGroupNames)
        groupCounts(g,k) = nnz(comparisonMask & R{k}.bestGroupIdx==g);
    end
    winnerPct(:,k) = winnerCounts(:,k)./nPixel.*100;
    categoryPct(:,k) = categoryCounts(:,k)./nPixel.*100;
    groupPct(:,k) = groupCounts(:,k)./nPixel.*100;
    if sum(winnerCounts(:,k)) ~= nPixel || ...
            sum(categoryCounts(:,k)) ~= nPixel || ...
            sum(groupCounts(:,k)) ~= nPixel
        error('Composition does not sum to shared domain for %s.', ...
            schemeNames{k});
    end
end

winnerComposition = table(modelNames,winnerCounts(:,1),winnerPct(:,1), ...
    winnerCounts(:,2),winnerPct(:,2),winnerCounts(:,3),winnerPct(:,3), ...
    winnerPct(:,2)-winnerPct(:,1),winnerPct(:,3)-winnerPct(:,1), ...
    'VariableNames',{'Model','M0_Count','M0_Pct','A_Count','A_Pct', ...
    'B_Count','B_Pct','A_minus_M0_pp','B_minus_M0_pp'});

categoryComposition = table(categoryNames,categoryCounts(:,1),categoryPct(:,1), ...
    categoryCounts(:,2),categoryPct(:,2),categoryCounts(:,3),categoryPct(:,3), ...
    categoryPct(:,2)-categoryPct(:,1),categoryPct(:,3)-categoryPct(:,1), ...
    'VariableNames',{'Category','M0_Count','M0_Pct','A_Count','A_Pct', ...
    'B_Count','B_Pct','A_minus_M0_pp','B_minus_M0_pp'});

groupComposition = table(modelGroupNames(:),groupCounts(:,1),groupPct(:,1), ...
    groupCounts(:,2),groupPct(:,2),groupCounts(:,3),groupPct(:,3), ...
    groupPct(:,2)-groupPct(:,1),groupPct(:,3)-groupPct(:,1), ...
    'VariableNames',{'Group','M0_Count','M0_Pct','A_Count','A_Pct', ...
    'B_Count','B_Pct','A_minus_M0_pp','B_minus_M0_pp'});

%% ================ 5. M0 consistency against original s2 ==============
% First, ALWAYS compare manuscript-rounded M0 category percentages.
% If historical s2 MAT files exist, also compare maps pixel by pixel.
expectedDiffPP = abs(categoryPct(:,1)-expectedM0ClassPct);
old5AuditAvailable = false;
old6AuditAvailable = false;
old5MaskMismatch = NaN;
old5WinnerMismatch = NaN;
old5CategoryMismatch = NaN;
old6CategoryMismatch = NaN;

if isfile(oldFiveClassFile)
    try
        SOld = load(oldFiveClassFile,'R_original');
        if ~isfield(SOld,'R_original')
            error('R_original variable was not found.');
        end
        oldR = SOld.R_original;
        fields = {'analysisMask','bestAll9Idx','gainClassMap'};
        for t = 1:numel(fields)
            if ~isfield(oldR,fields{t})
                error('R_original has no field %s.',fields{t});
            end
        end
        if ~isequal(size(oldR.gainClassMap),[nLat,nLon])
            error('Historical 5-class map has incompatible dimensions.');
        end
        old5AuditAvailable = true;
        old5MaskMismatch = nnz(xor(logical(oldR.analysisMask),comparisonMask));
        old5WinnerMismatch = nnz(comparisonMask & ...
            oldR.bestAll9Idx~=R_M0.bestAll9Idx);
        old5CategoryMismatch = nnz(comparisonMask & ...
            oldR.gainClassMap~=R_M0.gainClassMap);
    catch ME
        warning('Cannot audit original 5-class MAT: %s',ME.message);
    end
else
    fprintf('Original s2_gain_attribution2 MAT not found: optional audit skipped.\n');
end

if isfile(oldSixClassFile)
    try
        SOld6 = load(oldSixClassFile,'A_original');
        if ~isfield(SOld6,'A_original') || ...
                ~isfield(SOld6.A_original,'classMasks')
            error('Historical A_original.classMasks is unavailable.');
        end
        masksOld = SOld6.A_original.classMasks;
        if numel(masksOld) ~= nClass
            error('Historical classMasks length is not six.');
        end
        oldClassMap = nan(nLat,nLon);
        for c = 1:nClass
            oldClassMap(logical(masksOld{c})) = c;
        end
        old6CategoryMismatch = nnz(comparisonMask & ...
            oldClassMap~=R_M0.gainClass6Map);
        old6AuditAvailable = true;
    catch ME
        warning('Cannot audit original 6-class MAT: %s',ME.message);
    end
else
    fprintf('Original s2_gain_attribution3 MAT not found: optional audit skipped.\n');
end

M0AuditNames = { ...
    'Stage1_all_nine_models_PASS'; ...
    'M0_class_pct_matches_manuscript_1decimal'; ...
    'Old_5class_results_available'; ...
    'Old_5class_mask_mismatches'; ...
    'Old_5class_winner_mismatches'; ...
    'Old_5class_category_mismatches'; ...
    'Old_6class_results_available'; ...
    'Old_6class_category_mismatches'};
M0AuditValues = [1;double(all(expectedDiffPP<=expectedM0ClassTolerancePP)); ...
    double(old5AuditAvailable);old5MaskMismatch;old5WinnerMismatch; ...
    old5CategoryMismatch;double(old6AuditAvailable);old6CategoryMismatch];
M0Audit = table(M0AuditNames,M0AuditValues, ...
    'VariableNames',{'Check','Value'});

if any(expectedDiffPP>expectedM0ClassTolerancePP)
    error(['M0 detailed class percentages do not reproduce the current ', ...
        'manuscript rounded percentages. Check original rules.']);
end
if old5AuditAvailable && (old5MaskMismatch~=0 || ...
        old5WinnerMismatch~=0 || old5CategoryMismatch~=0)
    error('M0 differs from historical s2_gain_attribution2 results.');
end
if old6AuditAvailable && old6CategoryMismatch~=0
    error('M0 differs from historical s2_gain_attribution3 six-class maps.');
end
fprintf('M0 six-class percentages match manuscript rounding.\n');
if old5AuditAvailable
    fprintf('M0 vs original 5-class audit: mask=%d winner=%d class=%d mismatch(es).\n', ...
        old5MaskMismatch,old5WinnerMismatch,old5CategoryMismatch);
end
if old6AuditAvailable
    fprintf('M0 vs original 6-class audit: class=%d mismatch(es).\n', ...
        old6CategoryMismatch);
end
fprintf('\n');

%% ==================== 6. Pairwise agreement ==========================
winnerTransitionCounts = zeros(nModel,nModel,nPair);
categoryTransitionCounts = zeros(nClass,nClass,nPair);
winnerTransitionRowPct = nan(nModel,nModel,nPair);
categoryTransitionRowPct = nan(nClass,nClass,nPair);

winnerSameCount = zeros(nPair,1);
categorySameCount = zeros(nPair,1);
groupSameCount = zeros(nPair,1);

for p = 1:nPair
    ia = pairSchemeIdx(p,1);
    ib = pairSchemeIdx(p,2);
    baseWinner = R{ia}.bestAll9Idx(comparisonMask);
    altWinner = R{ib}.bestAll9Idx(comparisonMask);
    baseCat = R{ia}.gainClass6Map(comparisonMask);
    altCat = R{ib}.gainClass6Map(comparisonMask);
    baseGroup = R{ia}.bestGroupIdx(comparisonMask);
    altGroup = R{ib}.bestGroupIdx(comparisonMask);

    winnerSameCount(p) = nnz(baseWinner==altWinner);
    categorySameCount(p) = nnz(baseCat==altCat);
    groupSameCount(p) = nnz(baseGroup==altGroup);

    % Rows = reference scheme; columns = alternative scheme.
    for i = 1:nModel
        for j = 1:nModel
            winnerTransitionCounts(i,j,p) = nnz(baseWinner==i & altWinner==j);
        end
    end
    for i = 1:nClass
        for j = 1:nClass
            categoryTransitionCounts(i,j,p) = nnz(baseCat==i & altCat==j);
        end
    end
    for i = 1:nModel
        rowN = sum(winnerTransitionCounts(i,:,p));
        if rowN>0
            winnerTransitionRowPct(i,:,p) = ...
                winnerTransitionCounts(i,:,p)./rowN.*100;
        end
    end
    for i = 1:nClass
        rowN = sum(categoryTransitionCounts(i,:,p));
        if rowN>0
            categoryTransitionRowPct(i,:,p) = ...
                categoryTransitionCounts(i,:,p)./rowN.*100;
        end
    end
    if sum(sum(winnerTransitionCounts(:,:,p)))~=nPixel || ...
       sum(sum(categoryTransitionCounts(:,:,p)))~=nPixel
        error('Transition matrices have an invalid cell total for %s.', ...
            pairNames{p});
    end
end

agreement = table(pairNames(:),repmat(nPixel,nPair,1), ...
    winnerSameCount,nPixel-winnerSameCount, ...
    winnerSameCount./nPixel.*100, ...
    categorySameCount,nPixel-categorySameCount, ...
    categorySameCount./nPixel.*100, ...
    groupSameCount,nPixel-groupSameCount,groupSameCount./nPixel.*100, ...
    'VariableNames',{'Comparison','ValidPixels','WinnerSame', ...
    'WinnerChanged','WinnerAgreementPct', ...
    'CategorySame','CategoryChanged','CategoryAgreementPct', ...
    'GroupSame','GroupChanged','GroupAgreementPct'});

%% ======================== 7. Diagnostics ===========================
numExactWinnerTies = zeros(nScheme,1);
numExactGroupTies = zeros(nScheme,1);
numAssigned = zeros(nScheme,1);
numUnassigned = zeros(nScheme,1);
numMultiAssigned = zeros(nScheme,1);
numWinnerMismatch = zeros(nScheme,1);

for k = 1:nScheme
    numExactWinnerTies(k) = R{k}.nExactWinnerTies;
    numExactGroupTies(k) = R{k}.nExactGroupTies;
    numAssigned(k) = R{k}.nAssigned;
    numUnassigned(k) = R{k}.nUnassigned;
    numMultiAssigned(k) = R{k}.nMultiAssigned;
    numWinnerMismatch(k) = R{k}.nWinnerMismatch;
end

diagnostics = table(schemeNames,numAssigned,numUnassigned, ...
    numMultiAssigned,numWinnerMismatch,numExactWinnerTies, ...
    numExactGroupTies, ...
    'VariableNames',{'Scheme','Assigned','Unassigned','MultiAssigned', ...
    'WinnerReconstructionMismatches','ExactModelWinnerTies','ExactGroupTies'});

fprintf('======================= WINNER COMPOSITION =======================\n');
disp(winnerComposition);
fprintf('====================== CATEGORY COMPOSITION ======================\n');
disp(categoryComposition);
fprintf('======================== GROUP COMPOSITION =======================\n');
disp(groupComposition);
fprintf('========================== AGREEMENT ============================\n');
disp(agreement);
fprintf('========================== DIAGNOSTICS ===========================\n');
disp(diagnostics);

%% ========================== 8. Save ===============================
% The three R_* structs contain all the group fusion references and masks
% needed for the Stage 3 quantitative gain decomposition.
resultMatFile = fullfile(outDir,'s5_NRMSE_sensitivity2_results.mat');
AR_NRMSE_M0 = scoreMaps{1}(:,:,1);
AR_NRMSE_A = scoreMaps{2}(:,:,1);
AR_NRMSE_B = scoreMaps{3}(:,:,1);

save(resultMatFile, ...
    'R_M0','R_A','R_B', ...
    'AR_NRMSE_M0','AR_NRMSE_A','AR_NRMSE_B', ...
    'modelNames','schemeNames','categoryNames','modelGroupNames', ...
    'testYears','lat','lon','comparisonMask','M0ReferenceMask', ...
    'winnerCounts','winnerPct','categoryCounts','categoryPct', ...
    'groupCounts','groupPct','pairNames','pairSchemeIdx', ...
    'winnerTransitionCounts','winnerTransitionRowPct', ...
    'categoryTransitionCounts','categoryTransitionRowPct', ...
    'winnerComposition','categoryComposition','groupComposition', ...
    'agreement','diagnostics','M0Audit', ...
    'stage1File','oldFiveClassFile','oldSixClassFile', ...
    '-v7.3');

writetable(winnerComposition, ...
    fullfile(outDir,'s5_stage2_winner_composition.csv'));
writetable(categoryComposition, ...
    fullfile(outDir,'s5_stage2_category_composition.csv'));
writetable(groupComposition, ...
    fullfile(outDir,'s5_stage2_group_composition.csv'));
writetable(agreement,fullfile(outDir,'s5_stage2_agreement.csv'));
writetable(diagnostics,fullfile(outDir,'s5_stage2_diagnostics.csv'));
writetable(M0Audit,fullfile(outDir,'s5_stage2_M0_audit.csv'));

excelFile = fullfile(outDir,'s5_NRMSE_sensitivity2_tables.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(winnerComposition,excelFile,'Sheet','Winners');
writetable(categoryComposition,excelFile,'Sheet','Categories');
writetable(groupComposition,excelFile,'Sheet','Groups');
writetable(agreement,excelFile,'Sheet','Agreement');
writetable(diagnostics,excelFile,'Sheet','Diagnostics');
writetable(M0Audit,excelFile,'Sheet','M0_Audit');

for p = 1:nPair
    shortName = strrep(pairNames{p},'_vs_','_');
    wt = transitionTable(winnerTransitionCounts(:,:,p), ...
        modelNames,modelNames,'FromModel','To_');
    wpt = transitionTable(winnerTransitionRowPct(:,:,p), ...
        modelNames,modelNames,'FromModel','To_');
    ct = transitionTable(categoryTransitionCounts(:,:,p), ...
        categoryNames,{'NoGain';'Linear';'Nonlinear';'Hydro'; ...
        'MixedL';'MixedN'},'FromCategory','To_');
    cpt = transitionTable(categoryTransitionRowPct(:,:,p), ...
        categoryNames,{'NoGain';'Linear';'Nonlinear';'Hydro'; ...
        'MixedL';'MixedN'},'FromCategory','To_');

    prefix = ['s5_stage2_transitions_',shortName];
    writetable(wt,fullfile(outDir,[prefix,'_winner_counts.csv']));
    writetable(wpt,fullfile(outDir,[prefix,'_winner_rowpct.csv']));
    writetable(ct,fullfile(outDir,[prefix,'_category_counts.csv']));
    writetable(cpt,fullfile(outDir,[prefix,'_category_rowpct.csv']));

    writetable(wt,excelFile,'Sheet',['W_',shortName,'_Cnt']);
    writetable(wpt,excelFile,'Sheet',['W_',shortName,'_Pct']);
    writetable(ct,excelFile,'Sheet',['C_',shortName,'_Cnt']);
    writetable(cpt,excelFile,'Sheet',['C_',shortName,'_Pct']);
end

fprintf('\nSaved MAT: %s\n',resultMatFile);
fprintf('Saved Excel: %s\n',excelFile);
fprintf('Elapsed: %.1f s\n',toc(startTic));
fprintf('STAGE 2 COMPLETED: M0/A/B winners and six-class attribution maps.\n');
fprintf('NOTE: component contribution shares belong to Stage 3.\n');
fprintf('=============================================================\n');
diary off;

%% ======================== Local functions =============================
function R = classifyThreeYearWinners(modelStack,analysisMask)
% This is a literal adaptation of the original s2_gain_attribution2.m and
% the category split in s2_gain_attribution3.m. Score maps are 3-year mean
% NRMSE for the ORIGINAL run of nine models in the original model order.

    [nr,nc,nModel] = size(modelStack);
    assert(nModel==9,'Expected nine forecasting configurations.');
    assert(isequal(size(analysisMask),[nr,nc]), ...
        'Analysis mask and modelStack spatial size differ.');

    [directAll9NRMSE,directBestAll9Idx] = min(modelStack,[],3);

    linearStack = modelStack(:,:,1:3);
    nonlinearStack = modelStack(:,:,4:6);
    nonlinearXStack = modelStack(:,:,7:9);

    [linearFusion,bestLinearLocalIdx] = min(linearStack,[],3);
    [nonlinearFusion,bestNonlinearLocalIdx] = min(nonlinearStack,[],3);
    [nonlinearXFusion,bestNonlinearXLocalIdx] = min(nonlinearXStack,[],3);

    groupFusionStack = cat(3,linearFusion,nonlinearFusion,nonlinearXFusion);
    [sortedGroupNRMSE,sortedGroupIdx] = sort(groupFusionStack,3,'ascend');
    bestGroupIdx = sortedGroupIdx(:,:,1);
    secondGroupIdx = sortedGroupIdx(:,:,2);
    bestGroupNRMSE = sortedGroupNRMSE(:,:,1);

    % This grouping step reproduces the original nine-model selection.
    bestLinearGlobalIdx = bestLinearLocalIdx;
    bestNonlinearGlobalIdx = bestNonlinearLocalIdx+3;
    bestNonlinearXGlobalIdx = bestNonlinearXLocalIdx+6;
    bestAll9Idx = nan(nr,nc);
    maskBestLinear = analysisMask & bestGroupIdx==1;
    maskBestNonlinear = analysisMask & bestGroupIdx==2;
    maskBestNonlinearX = analysisMask & bestGroupIdx==3;
    bestAll9Idx(maskBestLinear) = bestLinearGlobalIdx(maskBestLinear);
    bestAll9Idx(maskBestNonlinear) = bestNonlinearGlobalIdx(maskBestNonlinear);
    bestAll9Idx(maskBestNonlinearX) = ...
        bestNonlinearXGlobalIdx(maskBestNonlinearX);

    % Paired antecedent-hydroclimatic reference must be the SAME
    % algorithm type: SVM-SVMX, RF-RFX, LSTM-LSTMX.
    maskHydroClean = analysisMask & bestGroupIdx==3 & ...
        secondGroupIdx==2 & ...
        bestNonlinearXLocalIdx==bestNonlinearLocalIdx;

    maskNoGain = analysisMask & bestGroupIdx==1 & bestLinearLocalIdx==1;
    maskLinearGain = analysisMask & bestGroupIdx==1 & ...
        (bestLinearLocalIdx==2 | bestLinearLocalIdx==3);
    maskNonlinearGain = analysisMask & bestGroupIdx==2;
    maskMixedLinearSecond = analysisMask & bestGroupIdx==3 & ...
        secondGroupIdx==1;
    maskMixedNonlinearSecond = analysisMask & bestGroupIdx==3 & ...
        secondGroupIdx==2 & ...
        bestNonlinearXLocalIdx~=bestNonlinearLocalIdx;
    maskMixedGain = maskMixedLinearSecond | maskMixedNonlinearSecond;

    classMasks = {maskNoGain;maskLinearGain;maskNonlinearGain; ...
        maskHydroClean;maskMixedLinearSecond;maskMixedNonlinearSecond};
    classSum = zeros(nr,nc);
    gainClass6Map = nan(nr,nc);
    for c = 1:numel(classMasks)
        classSum = classSum + double(classMasks{c});
        gainClass6Map(classMasks{c}) = c;
    end
    gainClassMap = gainClass6Map;
    gainClassMap(gainClass6Map==6) = 5; % Original 5-class convention.

    nAssigned = nnz(analysisMask & classSum==1);
    nUnassigned = nnz(analysisMask & classSum==0);
    nMultiAssigned = nnz(analysisMask & classSum>1);
    nWinnerMismatch = nnz(analysisMask & ...
        (bestAll9Idx~=directBestAll9Idx | ...
         abs(bestGroupNRMSE-directAll9NRMSE)>1e-10));

    % Exact ties are diagnostic only, NEVER resolved using new tolerance.
    nExactWinnerTies = nnz(analysisMask & ...
        sum(modelStack==directAll9NRMSE,3)>1);
    nExactGroupTies = nnz(analysisMask & ...
        sum(groupFusionStack==bestGroupNRMSE,3)>1);

    % Exclude all cells outside shared comparison domain from maps.
    linearFusion(~analysisMask)=NaN;
    nonlinearFusion(~analysisMask)=NaN;
    nonlinearXFusion(~analysisMask)=NaN;
    bestLinearLocalIdx(~analysisMask)=NaN;
    bestNonlinearLocalIdx(~analysisMask)=NaN;
    bestNonlinearXLocalIdx(~analysisMask)=NaN;
    bestLinearGlobalIdx(~analysisMask)=NaN;
    bestNonlinearGlobalIdx(~analysisMask)=NaN;
    bestNonlinearXGlobalIdx(~analysisMask)=NaN;
    bestGroupIdx(~analysisMask)=NaN;
    secondGroupIdx(~analysisMask)=NaN;
    bestGroupNRMSE(~analysisMask)=NaN;
    directAll9NRMSE(~analysisMask)=NaN;
    directBestAll9Idx(~analysisMask)=NaN;

    R.analysisMask = analysisMask;
    R.linearFusion = linearFusion;
    R.nonlinearFusion = nonlinearFusion;
    R.nonlinearXFusion = nonlinearXFusion;
    R.bestLinearLocalIdx = bestLinearLocalIdx;
    R.bestNonlinearLocalIdx = bestNonlinearLocalIdx;
    R.bestNonlinearXLocalIdx = bestNonlinearXLocalIdx;
    R.bestLinearGlobalIdx = bestLinearGlobalIdx;
    R.bestNonlinearGlobalIdx = bestNonlinearGlobalIdx;
    R.bestNonlinearXGlobalIdx = bestNonlinearXGlobalIdx;
    R.bestGroupIdx = bestGroupIdx;
    R.secondGroupIdx = secondGroupIdx;
    R.bestAll9Idx = bestAll9Idx;
    R.allFusion = bestGroupNRMSE;
    R.maskNoGain = maskNoGain;
    R.maskLinearGain = maskLinearGain;
    R.maskNonlinearGain = maskNonlinearGain;
    R.maskHydroClean = maskHydroClean;
    R.maskMixedGain = maskMixedGain;
    R.maskMixedLinearSecond = maskMixedLinearSecond;
    R.maskMixedNonlinearSecond = maskMixedNonlinearSecond;
    R.gainClassMap = gainClassMap;
    R.gainClass6Map = gainClass6Map;
    R.classMasks = classMasks;
    R.nAssigned = nAssigned;
    R.nUnassigned = nUnassigned;
    R.nMultiAssigned = nMultiAssigned;
    R.nWinnerMismatch = nWinnerMismatch;
    R.nExactWinnerTies = nExactWinnerTies;
    R.nExactGroupTies = nExactGroupTies;
end

function T = transitionTable(data,rowLabels,columnLabels,rowHeading,columnPrefix)
% Make a simple CSV/Excel-friendly transition table, including row labels.
    data = double(data);
    assert(size(data,1)==numel(rowLabels) && ...
           size(data,2)==numel(columnLabels), ...
           'Transition table dimensions do not match the label counts.');
    columnNames = cell(numel(columnLabels),1);
    for j = 1:numel(columnLabels)
        columnNames{j} = [columnPrefix,columnLabels{j}];
    end
    T = array2table(data,'VariableNames',columnNames);
    T = addvars(T,rowLabels(:),'Before',1,'NewVariableNames',rowHeading);
end
