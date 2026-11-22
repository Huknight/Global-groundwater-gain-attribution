%% s2_gain_attribution3.m

clc;
clearvars;

%
% ==================== 1. Code Description ====================
% Purpose:
%   Reproduce the manuscript-style pixel-wise hierarchical attribution
%   profile and Table 1 for three stochastic model runs:
%       1. Original run
%       2. UC1 run
%       3. UC2 run
%
% Input:
%   s2_gain_attribution2_results.mat
%
% Figures:
%   Fig1: Pixel-wise hierarchical attribution profile - Original run.
%   Fig2: Pixel-wise hierarchical attribution profile - UC1 run.
%   Fig3: Pixel-wise hierarchical attribution profile - UC2 run.
%
% Tables:
%   Table1_Original:
%       Winner-conditioned gain decomposition for the original run.
%
%   Table1_UC1:
%       Winner-conditioned gain decomposition for UC1.
%
%   Table1_UC2:
%       Winner-conditioned gain decomposition for UC2.
%
% Positive-gain attribution classes:
%   1. Linear gain
%   2. Nonlinear gain
%   3. Hydroclimatic gain
%   4. Mixed gain with linear second-best
%   5. Mixed gain with nonlinear second-best
%
% Gain components:
%   1. Linear component
%   2. Nonlinear component
%   3. Hydroclimatic component
%   4. Mixed component
%
% Figure design:
%   - Pixels are first grouped by attribution class.
%   - Pixels within each class are sorted by relative total NRMSE gain.
%   - No-gain pixels are not plotted.
%   - Horizontal dashed lines indicate 10%, 20%, 30%, and 40%.
%   - Vertical dashed lines indicate the P90 and P95 positions within
%     each attribution class.
%   - Class annotations report:
%         n   = pixel number
%         Red = mean total NRMSE reduction
%         Pct = mean total NRMSE reduction relative to AR
%
% Table 1:
%   Contribution percentages represent the share of the summed total NRMSE
%   reduction explained by each component within a given attribution class.
%   Mean component reductions are reported in NRMSE units.
%
% Outputs:
%   1. Three PNG figures.
%   2. Three Table 1 CSV files.
%   3. One Excel workbook containing the three Table 1 sheets.
%   4. Three figure-annotation summary CSV files.
%   5. One MAT result file.
%
% Data source:
%   The attribution masks and component inputs are read from the revised
%   s2_gain_attribution2_results.mat produced by s2_gain_attribution2.m
%   on the revised master database, in which the 0.25-degree to 1-degree
%   upscaling has been corrected with spherical-area weights. The AR
%   reference NRMSE is read from the revised AR result file.
%
% Notes:
%   1. Mixed gain is separated into Mixed-L and Mixed-N in this script
%      because this distinction is required by the manuscript Figure 7
%      and Table 1.
%
%   2. The AR model is deterministic and is shared by all three runs.
%
%   3. A common y-axis limit is used for Fig1-Fig3 to support direct
%      comparison between the three runs.
%
% ==================== 2. Centralized Parameter Settings ====================

% ------------------------------ Paths ------------------------------
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';

dataRoot = fullfile(projectRoot, 's1_master_database');

linearDir = fullfile(dataRoot, '04_linear_models');

% Previous gain-attribution results from the revised database
inputDir = fullfile(projectRoot, ...
    's3_darft', ...
    's1_analysis', ...
    's2_gain_attribution2_results');

inputMatFile = fullfile(inputDir, ...
    's2_gain_attribution2_results.mat');

outDir = fullfile(projectRoot, ...
    's3_darft', ...
    's1_analysis', ...
    's2_gain_attribution3_results');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ------------------------ AR result file ------------------------
arFile = fullfile(linearDir, ...
    'AR_univariate_global_200401_202412.mat');

% ------------------- Attribution class definitions -------------------
classNames = { ...
    'No gain'; ...
    'Linear gain'; ...
    'Nonlinear gain'; ...
    'Hydroclimatic gain'; ...
    'Mixed gain with linear second-best'; ...
    'Mixed gain with nonlinear second-best'};

positiveClassNames = { ...
    'Linear gain'; ...
    'Nonlinear gain'; ...
    'Hydroclimatic gain'; ...
    'Mixed gain with linear second-best'; ...
    'Mixed gain with nonlinear second-best'};

positiveClassShortNames = { ...
    'Linear'; ...
    'Nonlinear'; ...
    'Hydro'; ...
    'Mixed-L'; ...
    'Mixed-N'};

% ------------------- Component definitions and colors -------------------
componentLabels = { ...
    'Linear component', ...
    'Nonlinear component', ...
    'Hydroclimatic component', ...
    'Mixed component'};

componentColors = [ ...
    0.85, 0.10, 0.10;   % Linear component: red
    0.00, 0.35, 0.90;   % Nonlinear component: blue
    0.00, 0.60, 0.20;   % Hydroclimatic component: green
    1.00, 0.80, 0.00];  % Mixed component: yellow

% ---------------------- Unified figure settings ----------------------
settings.figPosition = [100, 100, 2000, 700];

settings.fontName = 'Times New Roman';

settings.titleFontSize = 16;
settings.axisLabelFontSize = 13;
settings.tickFontSize = 13;
settings.annotationFontSize = 16;
settings.legendFontSize = 16;
settings.referenceFontSize = 13;

% Top annotation position (move the n/Red/Pct text downward)
settings.annotationYFraction = 0.93;

% Percentile-label settings for 90th and 95th vertical lines
settings.percentileLabelFontSize = 8;
settings.percentileLabelRotation = 90;
settings.percentileLabelYFraction = 0.50;

settings.axisLineWidth = 0.8;
settings.classSeparatorLineWidth = 0.8;
settings.percentileLineWidth90 = 0.8;
settings.percentileLineWidth95 = 1.0;
settings.referenceLineWidth = 0.7;

settings.areaFaceAlpha = 0.90;

settings.classSeparatorColor = [0.20, 0.20, 0.20];
settings.percentile90Color = [0.40, 0.40, 0.40];
settings.percentile95Color = [0.10, 0.10, 0.10];
settings.referenceLineColor = [0.55, 0.55, 0.55];

settings.horizontalReferenceLevels = [10, 20, 30, 40];

settings.outputResolution = 300;

%
% ==================== 3. Main Body ====================

% ------------------- Load attribution results -------------------

% Variables loaded from the revised attribution-2 results:
%   R_original, R_uc1, R_uc2, validMask

requiredVariables = { ...
    'R_original', ...
    'R_uc1', ...
    'R_uc2', ...
    'validMask'};

S = load(inputMatFile, requiredVariables{:});

R_originalInput = S.R_original;
R_uc1Input = S.R_uc1;
R_uc2Input = S.R_uc2;

validMask = logical(S.validMask);

clear S

% ------------------- Load AR NRMSE -------------------
D_AR = load(arFile, 'AR_NRMSE_mean');

AR = double(D_AR.AR_NRMSE_mean);

clear D_AR

% --------------- Reconstruct original-run attribution ---------------
fprintf('\nReconstructing original-run gain components...\n');

A_original = reconstructAttributionComponents( ...
    R_originalInput, AR, classNames, 'Original');

% ------------------- Reconstruct UC1 attribution -------------------
fprintf('Reconstructing UC1 gain components...\n');

A_uc1 = reconstructAttributionComponents( ...
    R_uc1Input, AR, classNames, 'UC1');

% ------------------- Reconstruct UC2 attribution -------------------
fprintf('Reconstructing UC2 gain components...\n');

A_uc2 = reconstructAttributionComponents( ...
    R_uc2Input, AR, classNames, 'UC2');

% ------------------- Build ordered plotting profiles -------------------
P_original = buildAttributionProfile( ...
    A_original, positiveClassNames, positiveClassShortNames);

P_uc1 = buildAttributionProfile( ...
    A_uc1, positiveClassNames, positiveClassShortNames);

P_uc2 = buildAttributionProfile( ...
    A_uc2, positiveClassNames, positiveClassShortNames);

% ------------------- Determine common y-axis limit -------------------
maximumProfileGain = max([ ...
    P_original.yMaxData, ...
    P_uc1.yMaxData, ...
    P_uc2.yMaxData]);

% Round the maximum to the next 10-percentage-point level.
commonYBase = ceil(maximumProfileGain / 10) * 10;

% Ensure that the 40% reference line can always be displayed.
commonYBase = max(commonYBase, 40);

% Add upper space for class annotations.
commonYTop = commonYBase * 1.10;

% ------------------- Fig1: Original-run attribution profile -------------------
fig1File = fullfile(outDir, 'Fig1_Pixelwise_Attribution_Profile_Original.png');

plotAttributionProfile( ...
    P_original, componentLabels, componentColors, ...
    'Pixel-wise hierarchical attribution profile - Original run', ...
    commonYTop, fig1File, settings);

% ----------------------- Fig2: UC1 attribution profile -----------------------
fig2File = fullfile(outDir, 'Fig2_Pixelwise_Attribution_Profile_UC1.png');

plotAttributionProfile( ...
    P_uc1, componentLabels, componentColors, ...
    'Pixel-wise hierarchical attribution profile - UC1 run', ...
    commonYTop, fig2File, settings);

% ----------------------- Fig3: UC2 attribution profile -----------------------
fig3File = fullfile(outDir, 'Fig3_Pixelwise_Attribution_Profile_UC2.png');

plotAttributionProfile( ...
    P_uc2, componentLabels, componentColors, ...
    'Pixel-wise hierarchical attribution profile - UC2 run', ...
    commonYTop, fig3File, settings);

% ------------------- Table 1 outputs -------------------
Table1_Original = A_original.Table1;
Table1_UC1 = A_uc1.Table1;
Table1_UC2 = A_uc2.Table1;

table1OriginalCsv = fullfile(outDir, 'Table1_Attribution_Decomposition_Original.csv');
table1UC1Csv = fullfile(outDir, 'Table1_Attribution_Decomposition_UC1.csv');
table1UC2Csv = fullfile(outDir, 'Table1_Attribution_Decomposition_UC2.csv');

writetable(Table1_Original, table1OriginalCsv);
writetable(Table1_UC1, table1UC1Csv);
writetable(Table1_UC2, table1UC2Csv);

% ------------------- Figure-annotation summary tables -------------------
FigureSummary_Original = A_original.FigureSummaryTable;
FigureSummary_UC1 = A_uc1.FigureSummaryTable;
FigureSummary_UC2 = A_uc2.FigureSummaryTable;

figureSummaryOriginalCsv = fullfile(outDir, 'Figure7_Annotation_Summary_Original.csv');
figureSummaryUC1Csv = fullfile(outDir, 'Figure7_Annotation_Summary_UC1.csv');
figureSummaryUC2Csv = fullfile(outDir, 'Figure7_Annotation_Summary_UC2.csv');

writetable(FigureSummary_Original, figureSummaryOriginalCsv);
writetable(FigureSummary_UC1, figureSummaryUC1Csv);
writetable(FigureSummary_UC2, figureSummaryUC2Csv);

% ------------------- Excel workbook -------------------
excelFile = fullfile(outDir, 's2_gain_attribution3_tables.xlsx');

writetable(Table1_Original, excelFile, 'Sheet', 'Table1_Original');
writetable(Table1_UC1, excelFile, 'Sheet', 'Table1_UC1');
writetable(Table1_UC2, excelFile, 'Sheet', 'Table1_UC2');
writetable(FigureSummary_Original, excelFile, 'Sheet', 'FigSummary_Original');
writetable(FigureSummary_UC1, excelFile, 'Sheet', 'FigSummary_UC1');
writetable(FigureSummary_UC2, excelFile, 'Sheet', 'FigSummary_UC2');

% ------------------- Command-window output -------------------
displayRunResults('ORIGINAL RUN', A_original, P_original);
displayRunResults('UC1 RUN', A_uc1, P_uc1);
displayRunResults('UC2 RUN', A_uc2, P_uc2);

fprintf('\nCommon y-axis upper limit for Fig1-Fig3: %.2f%%\n', commonYTop);

% ------------------- Save MAT results -------------------
resultMatFile = fullfile(outDir, 's2_gain_attribution3_results.mat');

save(resultMatFile, ...
    'A_original', ...
    'A_uc1', ...
    'A_uc2', ...
    'P_original', ...
    'P_uc1', ...
    'P_uc2', ...
    'Table1_Original', ...
    'Table1_UC1', ...
    'Table1_UC2', ...
    'FigureSummary_Original', ...
    'FigureSummary_UC1', ...
    'FigureSummary_UC2', ...
    'classNames', ...
    'positiveClassNames', ...
    'positiveClassShortNames', ...
    'componentLabels', ...
    'componentColors', ...
    'commonYTop', ...
    'validMask', ...
    'settings', ...
    '-v7.3');

fprintf('\nCompleted: three-run hierarchical attribution profiles and Table 1.\n');
fprintf('Output directory:\n%s\n', outDir);

%
% ==================== 4. Additional Functions ====================

function A = reconstructAttributionComponents( ...
    R, AR, classNames, runName)
% Reconstruct the hierarchical gain components for one run from the
% attribution masks and the group-level fusion maps, and build the
% manuscript-style Table 1.

    %% ---------- Load basic arrays ----------

    analysisMask = logical(R.analysisMask);

    L = double(R.linearFusion);
    N = double(R.nonlinearFusion);
    X = double(R.nonlinearXFusion);
    B = double(R.allFusion);

    bestNonlinearLocalIdx = double(R.bestNonlinearLocalIdx);
    bestNonlinearXLocalIdx = double(R.bestNonlinearXLocalIdx);

    bestGroupIdx = double(R.bestGroupIdx);
    secondGroupIdx = double(R.secondGroupIdx);

    %% ---------- Attribution masks ----------

    maskNoGain = analysisMask & logical(R.maskNoGain);
    maskLinearGain = analysisMask & logical(R.maskLinearGain);
    maskNonlinearGain = analysisMask & logical(R.maskNonlinearGain);
    maskHydroClean = analysisMask & logical(R.maskHydroClean);

    % NonlinearX winner with Linear as the second-best group.
    maskMixedLinearSecond = ...
        analysisMask & ...
        bestGroupIdx == 3 & ...
        secondGroupIdx == 1;

    % NonlinearX winner with Nonlinear as the second-best group,
    % but the two algorithms are not paired.
    maskMixedNonlinearSecond = ...
        analysisMask & ...
        bestGroupIdx == 3 & ...
        secondGroupIdx == 2 & ...
        bestNonlinearXLocalIdx ~= bestNonlinearLocalIdx;

    classMasks = { ...
        maskNoGain; ...
        maskLinearGain; ...
        maskNonlinearGain; ...
        maskHydroClean; ...
        maskMixedLinearSecond; ...
        maskMixedNonlinearSecond};

    %% ---------- Class-assignment diagnostics ----------

    classSum = ...
        double(maskNoGain) + ...
        double(maskLinearGain) + ...
        double(maskNonlinearGain) + ...
        double(maskHydroClean) + ...
        double(maskMixedLinearSecond) + ...
        double(maskMixedNonlinearSecond);

    nValid = nnz(analysisMask);

    nAssigned = nnz(analysisMask & classSum == 1);
    nUnassigned = nnz(analysisMask & classSum == 0);
    nMultiAssigned = nnz(analysisMask & classSum > 1);

    %% ---------- Total gain ----------

    G_total = AR - B;
    G_total(~analysisMask) = NaN;

    %% ---------- Initialize gain components ----------

    G_linear = nan(size(AR));
    G_nonlinear = nan(size(AR));
    G_hydro = nan(size(AR));
    G_mixed = nan(size(AR));

    G_linear(analysisMask) = 0;
    G_nonlinear(analysisMask) = 0;
    G_hydro(analysisMask) = 0;
    G_mixed(analysisMask) = 0;

    %% ---------- Class 1: No gain ----------

    G_linear(maskNoGain) = 0;
    G_nonlinear(maskNoGain) = 0;
    G_hydro(maskNoGain) = 0;
    G_mixed(maskNoGain) = 0;

    %% ---------- Class 2: Linear gain ----------

    G_linear(maskLinearGain) = AR(maskLinearGain) - B(maskLinearGain);

    %% ---------- Class 3: Nonlinear gain ----------

    G_linear(maskNonlinearGain) = AR(maskNonlinearGain) - L(maskNonlinearGain);
    G_nonlinear(maskNonlinearGain) = L(maskNonlinearGain) - N(maskNonlinearGain);

    %% ---------- Class 4: Hydroclimatic gain ----------

    G_linear(maskHydroClean) = AR(maskHydroClean) - L(maskHydroClean);
    G_nonlinear(maskHydroClean) = L(maskHydroClean) - N(maskHydroClean);
    G_hydro(maskHydroClean) = N(maskHydroClean) - X(maskHydroClean);

    %% ---------- Class 5: Mixed-L ----------

    G_linear(maskMixedLinearSecond) = AR(maskMixedLinearSecond) - L(maskMixedLinearSecond);
    G_mixed(maskMixedLinearSecond) = L(maskMixedLinearSecond) - X(maskMixedLinearSecond);

    %% ---------- Class 6: Mixed-N ----------

    G_linear(maskMixedNonlinearSecond) = AR(maskMixedNonlinearSecond) - L(maskMixedNonlinearSecond);
    G_nonlinear(maskMixedNonlinearSecond) = L(maskMixedNonlinearSecond) - N(maskMixedNonlinearSecond);
    G_mixed(maskMixedNonlinearSecond) = N(maskMixedNonlinearSecond) - X(maskMixedNonlinearSecond);

    %% ---------- Remove tiny numerical noise ----------

    tolerance = 1e-10;

    G_total(abs(G_total) < tolerance) = 0;
    G_linear(abs(G_linear) < tolerance) = 0;
    G_nonlinear(abs(G_nonlinear) < tolerance) = 0;
    G_hydro(abs(G_hydro) < tolerance) = 0;
    G_mixed(abs(G_mixed) < tolerance) = 0;

    %% ---------- Reconstruction check ----------

    G_reconstructed = G_linear + G_nonlinear + G_hydro + G_mixed;

    reconstructionError = G_reconstructed - G_total;
    reconstructionError(~analysisMask) = NaN;

    maxAbsReconstructionError = max( ...
        abs(reconstructionError(analysisMask)), ...
        [], 'omitnan');

    %% ---------- Table 1 calculations ----------

    nClass = numel(classNames);

    LinearContributionPct = nan(nClass, 1);
    NonlinearContributionPct = nan(nClass, 1);
    HydroContributionPct = nan(nClass, 1);
    MixedContributionPct = nan(nClass, 1);

    MeanLinearReduction = nan(nClass, 1);
    MeanNonlinearReduction = nan(nClass, 1);
    MeanHydroReduction = nan(nClass, 1);
    MeanMixedReduction = nan(nClass, 1);

    PixelCount = zeros(nClass, 1);
    PixelPct = zeros(nClass, 1);
    MeanTotalReduction = nan(nClass, 1);
    MeanTotalReductionPct = nan(nClass, 1);

    for iClass = 1:nClass
        maskClass = classMasks{iClass};

        PixelCount(iClass) = nnz(maskClass);

        if nValid > 0
            PixelPct(iClass) = PixelCount(iClass) / nValid * 100;
        else
            PixelPct(iClass) = NaN;
        end

        if PixelCount(iClass) == 0
            continue;
        end

        totalValues = G_total(maskClass);
        linearValues = G_linear(maskClass);
        nonlinearValues = G_nonlinear(maskClass);
        hydroValues = G_hydro(maskClass);
        mixedValues = G_mixed(maskClass);

        MeanTotalReduction(iClass) = mean(totalValues, 'omitnan');

        validRelativeMask = maskClass & isfinite(AR) & AR > 0;

        if nnz(validRelativeMask) > 0
            MeanTotalReductionPct(iClass) = mean( ...
                G_total(validRelativeMask) ./ AR(validRelativeMask) .* 100, ...
                'omitnan');
        end

        MeanLinearReduction(iClass) = mean(linearValues, 'omitnan');
        MeanNonlinearReduction(iClass) = mean(nonlinearValues, 'omitnan');
        MeanHydroReduction(iClass) = mean(hydroValues, 'omitnan');
        MeanMixedReduction(iClass) = mean(mixedValues, 'omitnan');

        totalSumForRatio = sum(totalValues, 'omitnan');

        if abs(totalSumForRatio) > 1e-12
            LinearContributionPct(iClass) = ...
                sum(linearValues, 'omitnan') ./ totalSumForRatio .* 100;

            NonlinearContributionPct(iClass) = ...
                sum(nonlinearValues, 'omitnan') ./ totalSumForRatio .* 100;

            HydroContributionPct(iClass) = ...
                sum(hydroValues, 'omitnan') ./ totalSumForRatio .* 100;

            MixedContributionPct(iClass) = ...
                sum(mixedValues, 'omitnan') ./ totalSumForRatio .* 100;
        end
    end

    %% ---------- Manuscript-style Table 1 ----------

    Table1 = table( ...
        classNames, ...
        LinearContributionPct, ...
        NonlinearContributionPct, ...
        HydroContributionPct, ...
        MixedContributionPct, ...
        MeanLinearReduction, ...
        MeanNonlinearReduction, ...
        MeanHydroReduction, ...
        MeanMixedReduction, ...
        'VariableNames', { ...
        'AttributionClass', ...
        'LinearContributionPct', ...
        'NonlinearContributionPct', ...
        'HydroContributionPct', ...
        'MixedContributionPct', ...
        'MeanLinearReduction', ...
        'MeanNonlinearReduction', ...
        'MeanHydroReduction', ...
        'MeanMixedReduction'});

    %% ---------- Figure annotation summary ----------

    FigureSummaryTable = table( ...
        classNames, ...
        PixelCount, ...
        PixelPct, ...
        MeanTotalReduction, ...
        MeanTotalReductionPct, ...
        'VariableNames', { ...
        'AttributionClass', ...
        'Pixels', ...
        'PixelPct', ...
        'MeanTotalNRMSEReduction', ...
        'MeanTotalReductionPct'});

    %% ---------- Save output structure ----------

    A.runName = runName;

    A.AR_NRMSE = AR;
    A.analysisMask = analysisMask;
    A.nValid = nValid;

    A.classNames = classNames;
    A.classMasks = classMasks;

    A.maskNoGain = maskNoGain;
    A.maskLinearGain = maskLinearGain;
    A.maskNonlinearGain = maskNonlinearGain;
    A.maskHydroClean = maskHydroClean;
    A.maskMixedLinearSecond = maskMixedLinearSecond;
    A.maskMixedNonlinearSecond = maskMixedNonlinearSecond;

    A.G_total = G_total;
    A.G_linear = G_linear;
    A.G_nonlinear = G_nonlinear;
    A.G_hydro = G_hydro;
    A.G_mixed = G_mixed;

    A.Table1 = Table1;
    A.FigureSummaryTable = FigureSummaryTable;

    A.nAssigned = nAssigned;
    A.nUnassigned = nUnassigned;
    A.nMultiAssigned = nMultiAssigned;
    A.maxAbsReconstructionError = maxAbsReconstructionError;
end

function P = buildAttributionProfile( ...
    A, positiveClassNames, positiveClassShortNames)
% Order the positive-gain pixels by class and relative total gain, and
% build the stacked component matrix for the attribution-profile figure.

    %% ---------- Relative gain percentages ----------

    % The AR denominator is stored in A.AR_NRMSE by
    % reconstructAttributionComponents; percentages are relative to AR.
    ARsafe = double(A.AR_NRMSE);

    ARsafe(~A.analysisMask) = NaN;
    ARsafe(ARsafe <= 0) = NaN;

    GtotalPct = 100 .* A.G_total ./ ARsafe;
    GlinearPct = 100 .* A.G_linear ./ ARsafe;
    GnonlinearPct = 100 .* A.G_nonlinear ./ ARsafe;
    GhydroPct = 100 .* A.G_hydro ./ ARsafe;
    GmixedPct = 100 .* A.G_mixed ./ ARsafe;

    %% ---------- Remove tiny numerical noise ----------

    tolerance = 1e-10;

    GtotalPct(abs(GtotalPct) < tolerance) = 0;
    GlinearPct(abs(GlinearPct) < tolerance) = 0;
    GnonlinearPct(abs(GnonlinearPct) < tolerance) = 0;
    GhydroPct(abs(GhydroPct) < tolerance) = 0;
    GmixedPct(abs(GmixedPct) < tolerance) = 0;

    %% ---------- Positive-gain class masks ----------

    positiveClassMasks = { ...
        A.maskLinearGain; ...
        A.maskNonlinearGain; ...
        A.maskHydroClean; ...
        A.maskMixedLinearSecond; ...
        A.maskMixedNonlinearSecond};

    nPositiveClass = numel(positiveClassMasks);

    %% ---------- Build ordered pixel sequence ----------

    orderedIdx = [];

    classCounts = zeros(nPositiveClass, 1);
    classPct = zeros(nPositiveClass, 1);

    classMeanTotalReduction = nan(nPositiveClass, 1);
    classMeanTotalReductionPct = nan(nPositiveClass, 1);

    for iClass = 1:nPositiveClass
        maskClass = positiveClassMasks{iClass};

        idxClass = find(maskClass);

        totalGainClass = GtotalPct(idxClass);

        [~, orderClass] = sort(totalGainClass, 'ascend');

        idxClass = idxClass(orderClass);

        orderedIdx = [orderedIdx; idxClass]; %#ok<AGROW>

        classCounts(iClass) = numel(idxClass);

        if A.nValid > 0
            classPct(iClass) = classCounts(iClass) ./ A.nValid .* 100;
        else
            classPct(iClass) = NaN;
        end

        if classCounts(iClass) > 0
            classMeanTotalReduction(iClass) = ...
                mean(A.G_total(maskClass), 'omitnan');

            classMeanTotalReductionPct(iClass) = ...
                mean(GtotalPct(maskClass), 'omitnan');
        end
    end

    nPlot = numel(orderedIdx);

    if nPlot == 0
        error('No positive-gain pixels were found for %s.', A.runName);
    end

    %% ---------- Construct stacked component matrix ----------

    Ypct = [ ...
        GlinearPct(orderedIdx), ...
        GnonlinearPct(orderedIdx), ...
        GhydroPct(orderedIdx), ...
        GmixedPct(orderedIdx)];

    Ypct(~isfinite(Ypct)) = 0;

    totalPctOrdered = GtotalPct(orderedIdx);
    totalPctOrdered(~isfinite(totalPctOrdered)) = 0;

    reconstructedPct = sum(Ypct, 2);

    reconstructionErrorPct = reconstructedPct - totalPctOrdered;

    maxAbsReconstructionErrorPct = max( ...
        abs(reconstructionErrorPct), [], 'omitnan');

    %% ---------- Class positions ----------

    classEnds = cumsum(classCounts);

    classStarts = [1; classEnds(1:end-1) + 1];

    classCenters = (classStarts + classEnds) ./ 2;

    nonEmptyClassMask = classCounts > 0;

    labelCenters = classCenters(nonEmptyClassMask);
    labelShortNames = positiveClassShortNames(nonEmptyClassMask);

    labelCounts = classCounts(nonEmptyClassMask);
    labelPct = classPct(nonEmptyClassMask);

    labelMeanReduction = classMeanTotalReduction(nonEmptyClassMask);
    labelMeanReductionPct = classMeanTotalReductionPct(nonEmptyClassMask);

    %% ---------- Class separators ----------

    separatorX = classEnds(1:end-1) + 0.5;

    separatorX = unique(separatorX);

    separatorX = separatorX(separatorX > 1 & separatorX < nPlot);

    %% ---------- P90 and P95 positions ----------

    percentile90X = nan(nPositiveClass, 1);
    percentile95X = nan(nPositiveClass, 1);

    for iClass = 1:nPositiveClass
        nClassPixel = classCounts(iClass);

        if nClassPixel <= 0
            continue;
        end

        localPosition90 = max(1, ceil(0.90 .* nClassPixel));
        localPosition95 = max(1, ceil(0.95 .* nClassPixel));

        percentile90X(iClass) = classStarts(iClass) + localPosition90 - 1;
        percentile95X(iClass) = classStarts(iClass) + localPosition95 - 1;
    end

    %% ---------- Maximum plotted gain ----------

    yMaxData = max(reconstructedPct, [], 'omitnan');

    if isempty(yMaxData) || ~isfinite(yMaxData) || yMaxData <= 0
        yMaxData = 1;
    end

    %% ---------- Save plotting structure ----------

    P.runName = A.runName;

    P.orderedIdx = orderedIdx;
    P.nPlot = nPlot;

    P.Ypct = Ypct;
    P.totalPctOrdered = totalPctOrdered;

    P.classNames = positiveClassNames;
    P.classShortNames = positiveClassShortNames;

    P.classCounts = classCounts;
    P.classPct = classPct;

    P.classStarts = classStarts;
    P.classEnds = classEnds;
    P.classCenters = classCenters;

    P.labelCenters = labelCenters;
    P.labelShortNames = labelShortNames;
    P.labelCounts = labelCounts;
    P.labelPct = labelPct;
    P.labelMeanReduction = labelMeanReduction;
    P.labelMeanReductionPct = labelMeanReductionPct;

    P.separatorX = separatorX;
    P.percentile90X = percentile90X;
    P.percentile95X = percentile95X;

    P.yMaxData = yMaxData;

    P.maxAbsReconstructionErrorPct = maxAbsReconstructionErrorPct;
end

function plotAttributionProfile( ...
    P, componentLabels, componentColors, ...
    figureTitle, yTop, figureFile, settings)
% Plot the stacked pixel-wise attribution profile and export it to a PNG
% file.

    %% ---------- Create figure ----------

    fig = figure('Color', 'w', 'Position', settings.figPosition);

    ax = axes('Parent', fig);

    x = 1:P.nPlot;

    %% ---------- Plot stacked attribution profile ----------

    areaHandles = area(ax, x, P.Ypct, 'LineStyle', 'none');

    for iComponent = 1:numel(areaHandles)
        areaHandles(iComponent).FaceColor = componentColors(iComponent, :);
        areaHandles(iComponent).EdgeColor = 'none';
        areaHandles(iComponent).FaceAlpha = settings.areaFaceAlpha;
    end

    hold(ax, 'on');

    %% ---------- Horizontal reference levels ----------

    for iLevel = 1:numel(settings.horizontalReferenceLevels)
        thisLevel = settings.horizontalReferenceLevels(iLevel);

        if thisLevel >= yTop
            continue;
        end

        plot(ax, [1, P.nPlot], [thisLevel, thisLevel], '--', ...
            'Color', settings.referenceLineColor, ...
            'LineWidth', settings.referenceLineWidth, ...
            'HandleVisibility', 'off');
    end

    %% ---------- Class separators ----------

    for iSeparator = 1:numel(P.separatorX)
        thisX = P.separatorX(iSeparator);

        plot(ax, [thisX, thisX], [0, yTop], '-', ...
            'Color', settings.classSeparatorColor, ...
            'LineWidth', settings.classSeparatorLineWidth, ...
            'HandleVisibility', 'off');
    end

    %% ---------- P90 positions ----------

    for iClass = 1:numel(P.percentile90X)
        thisX = P.percentile90X(iClass);

        if ~isfinite(thisX)
            continue;
        end

        plot(ax, [thisX, thisX], [0, yTop], '--', ...
            'Color', settings.percentile90Color, ...
            'LineWidth', settings.percentileLineWidth90, ...
            'HandleVisibility', 'off');

        text(ax, thisX, yTop * settings.percentileLabelYFraction, '90th', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'Rotation', settings.percentileLabelRotation, ...
            'FontName', settings.fontName, ...
            'FontSize', settings.percentileLabelFontSize, ...
            'Color', settings.percentile90Color);
    end

    %% ---------- P95 positions ----------

    for iClass = 1:numel(P.percentile95X)
        thisX = P.percentile95X(iClass);

        if ~isfinite(thisX)
            continue;
        end

        plot(ax, [thisX, thisX], [0, yTop], '--', ...
            'Color', settings.percentile95Color, ...
            'LineWidth', settings.percentileLineWidth95, ...
            'HandleVisibility', 'off');

        text(ax, thisX, yTop * settings.percentileLabelYFraction, '95th', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'Rotation', settings.percentileLabelRotation, ...
            'FontName', settings.fontName, ...
            'FontSize', settings.percentileLabelFontSize, ...
            'Color', settings.percentile95Color);
    end

    %% ---------- Class annotations ----------

    annotationY = yTop .* settings.annotationYFraction;

    for iLabel = 1:numel(P.labelCenters)
        annotationText = sprintf( ...
            'n=%d\nRed=%.4f\nPct=%.1f%%', ...
            P.labelCounts(iLabel), ...
            P.labelMeanReduction(iLabel), ...
            P.labelMeanReductionPct(iLabel));

        text(ax, P.labelCenters(iLabel), annotationY, annotationText, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontName', settings.fontName, ...
            'FontSize', settings.annotationFontSize, ...
            'Color', [0, 0, 0]);
    end

    hold(ax, 'off');

    %% ---------- Axis settings ----------

    if P.nPlot == 1
        xlim(ax, [0.5, 1.5]);
    else
        xlim(ax, [1, P.nPlot]);
    end

    ylim(ax, [0, yTop]);

    set(ax, ...
        'XTick', P.labelCenters, ...
        'XTickLabel', P.labelShortNames, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.tickFontSize, ...
        'LineWidth', settings.axisLineWidth);

    box(ax, 'on');

    xlabel(ax, 'Pixels ordered by attribution class and total NRMSE gain', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    ylabel(ax, 'NRMSE gain relative to AR (%)', ...
        'FontName', settings.fontName, ...
        'FontSize', settings.axisLabelFontSize);

    title(ax, figureTitle, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.titleFontSize, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none');

    %% ---------- Legend ----------

    legendHandle = legend( ...
        ax, areaHandles(:), componentLabels, ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal');

    set(legendHandle, ...
        'FontName', settings.fontName, ...
        'FontSize', settings.legendFontSize, ...
        'Box', 'off');

    %% ---------- Export figure ----------

    exportgraphics(fig, figureFile, 'Resolution', settings.outputResolution);
end

function displayRunResults(displayName, A, P)
% Print the attribution diagnostics, class counts, Table 1, and figure
% annotation summary to the command window.

    fprintf('\n');
    fprintf('============================================================================\n');
    fprintf('%s\n', displayName);
    fprintf('============================================================================\n');

    fprintf('Valid pixels                         : %d\n', A.nValid);
    fprintf('Assigned pixels                      : %d\n', A.nAssigned);
    fprintf('Unassigned pixels                    : %d\n', A.nUnassigned);
    fprintf('Multi-assigned pixels                : %d\n', A.nMultiAssigned);

    fprintf('No-gain pixels                       : %d (%.2f%%)\n', ...
        nnz(A.maskNoGain), nnz(A.maskNoGain) ./ A.nValid .* 100);

    fprintf('Positive-gain pixels plotted         : %d (%.2f%%)\n', ...
        P.nPlot, P.nPlot ./ A.nValid .* 100);

    fprintf('Maximum absolute reconstruction error: %.12e\n', ...
        A.maxAbsReconstructionError);

    fprintf('Maximum percentage reconstruction err: %.12e\n', ...
        P.maxAbsReconstructionErrorPct);

    fprintf('\nPositive-gain classes:\n');

    for iClass = 1:numel(P.classNames)
        fprintf('%-42s : %6d pixels (%6.2f%%)\n', ...
            P.classNames{iClass}, ...
            P.classCounts(iClass), ...
            P.classPct(iClass));
    end

    fprintf('\nTable 1:\n');
    disp(A.Table1);

    fprintf('\nFigure annotation summary:\n');
    disp(A.FigureSummaryTable);
end
