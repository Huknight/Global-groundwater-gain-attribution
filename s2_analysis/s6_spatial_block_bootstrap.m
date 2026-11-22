%% s6_spatial_block_bootstrap.m
% =========================================================================
% PURPOSE:
%   Paired spatial-block bootstrap for global groundwater forecasting.
%   Address small between-model NRMSE differences and spatial uncertainty
%   in the six winner-conditioned gain categories and gain components.
%
% IMPORTANT SCOPE:
%   - ORIGINAL nine-model run, ORIGINAL test-year-range NRMSE (M0) only.
%   - Each grid cell's 2022-2024 mean NRMSE is kept intact.
%   - Reuse stored forecasts/attribution outputs; NO model retraining.
%   - Spatial blocks are sampled with replacement as WHOLE clusters.
%   - Every bootstrap draw is shared by all nine models, all classes,
%     all components, and both global weighting definitions (PAIRED).
%   - Class labels and gain components are NOT refitted after resampling.
%   - Main published statistics remain GRID-EQUAL (unweighted).
%     AREA-WEIGHTED values are reported as auxiliary sensitivity results.
%
% SPATIAL SCALES:   2, 3, and 5 geographic degrees (fixed-origin blocks).
% REPLICATIONS:     2,000 per block scale (user-adjustable below).
% UNCERTAINTY:      Two-sided 95% percentile bootstrap intervals.
%                    RFX is paired against ALL eight alternative models.
%                    For these eight simultaneous comparisons, a two-sided
%                    Bonferroni-adjusted percentile interval is also given
%                    (familywise 95%; individual 99.375% for eight tests).
%                    A zero-exclusion flag is provided for each contrast;
%                    it is NOT a claim of temporal/general causal inference.
%
% INPUT FILES (existing output folders):
%   s5_NRMSE_sensitivity1_results/s5_NRMSE_sensitivity1_results.mat
%       modelNames, testYears, lat, lon, comparisonMask, NRMSE_mean_M0
%   s2_gain_attribution3_results/s2_gain_attribution3_results.mat
%       A_original (classMasks; G_linear; G_nonlinear; G_hydro;
%                   G_mixed; G_total; AR_NRMSE; analysisMask)
%   Optional audit only:
%   s2_gain_attribution2_results/s2_gain_attribution2_results.mat
%       lat, lon, validMask
%
% OUTPUTS:
%   s6_spatial_block_bootstrap_results.mat
%   s6_spatial_block_bootstrap_tables.xlsx
%   s6_global_model_means.csv
%   s6_global_category_shares.csv
%   s6_global_hydro_components.csv
%   s6_model_means_CI.csv
%   s6_paired_model_differences_CI.csv (all eight comparisons; 95% and
%       Bonferroni-adjusted intervals; three scales, two weightings)
%   s6_paired_model_differences_3deg.csv (main equal-grid comparisons)
%   s6_gain_categories_CI.csv
%   s6_hydro_components_CI.csv
%   s6_gain_contrasts_CI.csv
%   s6_spatial_block_diagnostics.csv
%   s6_spatial_block_bootstrap_log.txt
%   Fig1_Paired_Model_Differences_CI.png (if makeFigures = true)
%   Fig2_Attribution_Uncertainty_CI.png (if makeFigures = true)
%   Fig3_RFX_All_Model_Contrasts_3deg.png (if makeFigures = true)
%
% INTERPRETATION:
%   These intervals measure sensitivity to clustered geographic
%   resampling under the chosen fixed block sizes. They do NOT quantify
%   uncertainty across independent future years, or correct the
%   retrospective winner-selection bias. Geographic degree blocks are
%   not equal-area blocks, especially at high latitudes.
% =========================================================================

clc;
clearvars;

%% ======================= 1. Centralized settings ======================
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
baseDir = fullfile(projectRoot, 's3_darft', 's1_analysis');

stage1File = fullfile(baseDir, 's5_NRMSE_sensitivity1_results', ...
    's5_NRMSE_sensitivity1_results.mat');
attributionFile = fullfile(baseDir, 's2_gain_attribution3_results', ...
    's2_gain_attribution3_results.mat');
optionalCoordFile = fullfile(baseDir, 's2_gain_attribution2_results', ...
    's2_gain_attribution2_results.mat');

outDir = fullfile(baseDir, 's6_spatial_block_bootstrap_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Bootstrap settings (same 2000 repetitions at EACH scale)
blockDegrees = [2, 3, 5];
nBootstrap = 2000;
confidenceLevel = 95;
randomSeed = 20261116;
makeFigures = true;

% Reference-data safeguards
expectedValidCells = 17284;
expectedCategoryCounts = [1406; 4343; 4503; 2083; 3008; 1941];
referenceNRMSE = [ ...
    0.25798480; 0.26488363; 0.28895773; ...
    0.27803497; 0.26015954; 0.26360296; ...
    0.27951991; 0.25512264; 0.26843248];
referenceNRMSETolerance = 1e-6;
componentTolerance = 1e-8;

modelNamesExpected = { ...
    'AR'; 'ARIMA'; 'SARIMA'; ...
    'SVM'; 'RF'; 'LSTM'; ...
    'SVMX'; 'RFX'; 'LSTMX'};
categoryNames = { ...
    'No gain'; 'Linear gain'; 'Nonlinear gain'; ...
    'Hydroclimatic gain'; 'Mixed-L'; 'Mixed-N'};
componentNames = { ...
    'Linear component'; 'Nonlinear component'; ...
    'History (Linear + Nonlinear)'; 'Hydroclimatic component'};
contrastNames = { ...
    'Structural minus Hydroclimatic category (pp)'; ...
    'History minus Hydroclimatic contribution (pp)'};
weightingNames = {'GridEqual'; 'AreaWeighted'};

nModel = numel(modelNamesExpected);
% Preserve the two originally reported pairs first, then append the other
% six comparisons. Never independently resample the competing models.
rfxIndex = find(strcmp(modelNamesExpected,'RFX'));
pairModelIndices = [1,5,2,3,4,6,7,9];
assert(isscalar(rfxIndex) && rfxIndex == 8 && ...
    numel(pairModelIndices) == nModel-1 && ...
    isequal(sort(pairModelIndices),setdiff(1:nModel,rfxIndex)), ...
    'RFX comparison indices must contain all eight other models once.');
nPair = numel(pairModelIndices);
pairNames = cell(nPair,1);
for p = 1:nPair
    pairNames{p} = sprintf('RFX minus %s',modelNamesExpected{pairModelIndices(p)});
end
% Eight two-sided comparisons form ONE prespecified model family.
% Keep 95% intervals for each pair and additionally report simultaneous
% Bonferroni intervals at level 100*(1-0.05/8) = 99.375%.
bonferroniConfidenceLevel = 100*(1-(1-confidenceLevel/100)/nPair);
nClass = numel(categoryNames);
nComponent = numel(componentNames);
nWeighting = numel(weightingNames);

% Figure settings, consistent with earlier scripts
settings.figPosition = [100, 100, 1200, 520];
settings.fontName = 'Times New Roman';
settings.axisLabelFontSize = 12;
settings.tickFontSize = 11;
settings.legendFontSize = 10;
settings.outputResolution = 300;

logFile = fullfile(outDir, 's6_spatial_block_bootstrap_log.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
startTic = tic;

fprintf('==============================================================\n');
fprintf('S6: PAIRED SPATIAL BLOCK BOOTSTRAP (ORIGINAL M0)\n');
fprintf('==============================================================\n');
fprintf('Spatial block sizes : %s degrees\n', mat2str(blockDegrees));
fprintf('Bootstrap repetitions per scale: %d\n', nBootstrap);
fprintf('Confidence intervals: %.1f%% percentile\n', confidenceLevel);
fprintf('RFX versus other models: %d paired comparisons\n',nPair);
fprintf('Bonferroni simultaneous familywise level: %.1f%%\n',confidenceLevel);
fprintf('Bonferroni individual percentile level : %.3f%%\n', ...
    bonferroniConfidenceLevel);
fprintf('Random seed          : %d (+ block degrees)\n', randomSeed);
fprintf('Model retraining     : NONE\n\n');

assert(all(blockDegrees > 0 & blockDegrees == floor(blockDegrees)), ...
    'Block sizes must be positive integer degree values.');
assert(nBootstrap >= 100, 'At least 100 bootstrap repetitions are required.');
assert(confidenceLevel > 0 && confidenceLevel < 100, ...
    'Confidence level must lie between 0 and 100.');

%% ========================== 2. Load inputs ===========================
checkMATVariables(stage1File, { ...
    'modelNames','testYears','lat','lon', ...
    'comparisonMask','NRMSE_mean_M0'});
checkMATVariables(attributionFile, {'A_original'});

S = load(stage1File, 'modelNames','testYears','lat','lon', ...
    'comparisonMask','NRMSE_mean_M0');
Q = load(attributionFile, 'A_original');
A = Q.A_original;
clear Q;

assert(isequal(S.modelNames(:), modelNamesExpected(:)), ...
    'Stage-1 model order is inconsistent with the original manuscript.');
assert(isequal(double(S.testYears(:)'), [2022 2023 2024]), ...
    'Expected test years are [2022 2023 2024].');

lat = double(S.lat(:));
lon = double(S.lon(:));
validMask = logical(S.comparisonMask);
[nLat,nLon] = size(validMask);
assert(nLat == numel(lat) && nLon == numel(lon), ...
    'The Stage-1 latitude/longitude vectors do not match the mask.');
assert(isequal(size(S.NRMSE_mean_M0), [nLat nLon nModel]), ...
    'NRMSE_mean_M0 must have size [nLat,nLon,9].');

requiredAFields = { ...
    'analysisMask','AR_NRMSE','classMasks', ...
    'G_total','G_linear','G_nonlinear','G_hydro','G_mixed'};
for j = 1:numel(requiredAFields)
    assert(isfield(A,requiredAFields{j}), ...
        'A_original is missing field: %s', requiredAFields{j});
end
assert(isequal(size(A.analysisMask), size(validMask)), ...
    'Attribution map shape does not match Stage-1 mask.');
assert(isequal(logical(A.analysisMask), validMask), ...
    'The Stage-1 and attribution-3 valid-cell masks differ.');
assert(iscell(A.classMasks) && numel(A.classMasks) == nClass, ...
    'Expected exactly six attribution class masks.');

% Cross-check geographic coordinate ordering when original attribution-2
% results are available. This is an audit, not an additional dependency.
if isfile(optionalCoordFile)
    coordFields = who('-file', optionalCoordFile);
    if all(ismember({'lat','lon','validMask'}, coordFields))
        C0 = load(optionalCoordFile, 'lat','lon','validMask');
        assert(isequal(double(C0.lat(:)),lat) && ...
               isequal(double(C0.lon(:)),lon), ...
            'Stage-1 and attribution-2 geographic coordinates differ.');
        assert(isequal(logical(C0.validMask),validMask), ...
            'Stage-1 and attribution-2 masks differ.');
        fprintf('Optional geographic-coordinate audit: PASS\n');
        clear C0;
    else
        warning('Optional attribution-2 MAT lacks coordinate audit fields.');
    end
else
    fprintf('Optional attribution-2 coordinate file not found; skipping audit.\n');
end

%% ===================== 3. Assemble valid pixels ======================
validIndex = find(validMask);
nValid = numel(validIndex);
assert(nValid == expectedValidCells, ...
    'Expected %d valid cells, found %d. Verify the data revision.', ...
    expectedValidCells, nValid);

% MATLAB's column-major order is used consistently for all maps.
[Lon2D,Lat2D] = meshgrid(lon,lat);
latValid = Lat2D(validIndex);
lonValid = Lon2D(validIndex);
clear Lon2D Lat2D;

allModelMaps = reshape(double(S.NRMSE_mean_M0),[],nModel);
modelNRMSE = allModelMaps(validIndex,:);
assert(all(isfinite(modelNRMSE(:))), ...
    'Some valid cells have nonfinite model NRMSE scores.');

% Original model-performance reference audit (from s5 log).
originalModelMeans = mean(modelNRMSE,1)';
maxReferenceError = max(abs(originalModelMeans-referenceNRMSE));
assert(maxReferenceError < referenceNRMSETolerance, ...
    'M0 model means failed the reference audit (max error %.3g).', ...
    maxReferenceError);

% Six disjoint category indicators, without relabeling any grid cell.
categoryIndicator = zeros(nValid,nClass);
for c = 1:nClass
    classMap = A.classMasks{c};
    assert(isequal(size(classMap),[nLat nLon]), ...
        'Category %d map has incorrect dimensions.',c);
    categoryIndicator(:,c) = double(logical(classMap(validIndex)));
end
assert(all(sum(categoryIndicator,2)==1), ...
    'Every valid cell must belong to exactly one of six categories.');
actualCategoryCounts = sum(categoryIndicator,1)';
assert(isequal(actualCategoryCounts,expectedCategoryCounts), ...
    'M0 category counts do not match the manuscript: %s', ...
    mat2str(actualCategoryCounts'));

% Groundwater-history / exogenous contributions are only reported inside
% the Hydroclimatic-gain category (fourth category).
AR = double(A.AR_NRMSE(validIndex));
G_total = double(A.G_total(validIndex));
G_linear = double(A.G_linear(validIndex));
G_nonlinear = double(A.G_nonlinear(validIndex));
G_hydro = double(A.G_hydro(validIndex));
G_mixed = double(A.G_mixed(validIndex));
assert(all(isfinite([AR;G_total;G_linear;G_nonlinear;G_hydro;G_mixed])), ...
    'Nonfinite values found in the original gain decomposition.');

% Cross-check that the attribution reference is identical to model M0.
maxARDifference = max(abs(AR-modelNRMSE(:,1)));
minimum9 = min(modelNRMSE,[],2);
maxTotalGainDifference = max(abs(G_total-(AR-minimum9)));
maxReconstructionError = max(abs( ...
    G_total-G_linear-G_nonlinear-G_hydro-G_mixed));
assert(maxARDifference < componentTolerance, ...
    'AR NRMSE differs between Stage-1 and attribution-3 files.');
assert(maxTotalGainDifference < componentTolerance, ...
    'Winner gain differs between Stage-1 and attribution-3 files.');
assert(maxReconstructionError < componentTolerance, ...
    'Gain-component telescoping-sum audit failed.');

hydroClass = categoryIndicator(:,4);
gainHydroClass = [G_linear,G_nonlinear,G_hydro,G_total] ...
    .* repmat(hydroClass,1,4);
assert(sum(gainHydroClass(:,4)) > 0, ...
    'Hydroclimatic-gain class must have positive summed total gain.');

% Equal-area longitude width (1 degree) and latitude width (1 degree):
% cell area is exactly proportional to cos(latitude of its centre).
% The common multiplicative spherical-area constant cancels in means.
areaWeight = cosd(latValid);
assert(all(isfinite(areaWeight) & areaWeight > 0), ...
    'Nonpositive or nonfinite 1-degree spherical-area weights.');

% Build sufficient statistics for a block-resampled ratio estimate.
% Each weighting mode has 20 columns:
%   1: total weight
%   2-10: weighted model NRMSE sums (nine models)
%   11-16: weighted six-class counts
%   17-20: weighted Hydro-class linear, nonlinear, hydro, total gain sums
% Both modes use EXACTLY the same resampled blocks.
baseFeatures = [ones(nValid,1), modelNRMSE, ...
    categoryIndicator, gainHydroClass];
nFeature = size(baseFeatures,2);
assert(nFeature == 20, 'Unexpected number of sufficient statistics.');
allFeatures = [baseFeatures, bsxfun(@times,baseFeatures,areaWeight)];

pointResults(1) = extractStatistics( ...
    sum(baseFeatures,1),rfxIndex,pairModelIndices);
pointResults(2) = extractStatistics( ...
    sum(bsxfun(@times,baseFeatures,areaWeight),1), ...
    rfxIndex,pairModelIndices);

fprintf('Original valid cells        : %d\n',nValid);
fprintf('Category counts             : %s\n',mat2str(actualCategoryCounts'));
fprintf('NRMSE M0 mean audit error   : %.3g\n',maxReferenceError);
fprintf('AR M0 map audit error       : %.3g\n',maxARDifference);
fprintf('Total gain audit error      : %.3g\n',maxTotalGainDifference);
fprintf('Gain reconstruction error   : %.3g\n',maxReconstructionError);
fprintf('Original unweighted RFX-AR  : %.8f\n',pointResults(1).pairedDifferences(1));
fprintf('Original unweighted RFX-RF  : %.8f\n',pointResults(1).pairedDifferences(2));
fprintf('Structural-gain grids       : %.4f%%\n',pointResults(1).structuralPct);
fprintf('Hydroclimatic-gain grids    : %.4f%%\n',pointResults(1).categoryPct(4));
fprintf('History/Hydro components    : %.4f%% / %.4f%%\n\n', ...
    pointResults(1).componentPct(3),pointResults(1).componentPct(4));

%% ======================= 4. Global point tables =======================
modelTable = table( ...
    modelNamesExpected(:), ...
    pointResults(1).modelMeans(:), ...
    pointResults(2).modelMeans(:), ...
    'VariableNames',{'Model','GridEqualMeanNRMSE','AreaWeightedMeanNRMSE'});

categoryTable = table( ...
    categoryNames(:),actualCategoryCounts(:), ...
    pointResults(1).categoryPct(:), ...
    pointResults(2).categoryPct(:), ...
    'VariableNames',{'GainCategory','CellCount', ...
    'GridEqualPct','AreaWeightedPct'});

componentTable = table( ...
    componentNames(:), ...
    pointResults(1).componentPct(:), ...
    pointResults(2).componentPct(:), ...
    'VariableNames',{'HydroClassComponent', ...
    'GridEqualPct','AreaWeightedPct'});

%% ======================= 5. Spatial block bootstrap ==================
% Cell array rows are converted into CSV/Excel tables after all scales.
meanRows = cell(0,8);
pairRows = cell(0,15);
categoryRows = cell(0,8);
componentRows = cell(0,8);
contrastRows = cell(0,9);
blockRows = cell(0,8);

bootstrapResults = repmat(struct(),numel(blockDegrees),1);

for iScale = 1:numel(blockDegrees)
    blockDeg = blockDegrees(iScale);
    fprintf('----------------- BLOCK SIZE: %d degrees -----------------\n',blockDeg);

    [blockID,nBlocks] = buildSpatialBlocks(latValid,lonValid,blockDeg);
    blockSizes = accumarray(blockID,1,[nBlocks,1],@sum,0);

    % Sparse aggregation avoids repeatedly expanding individual cells.
    blockMembership = sparse(blockID,(1:nValid)',1,nBlocks,nValid);
    blockFeatures = full(blockMembership * allFeatures);
    clear blockMembership;
    assert(max(abs(sum(blockFeatures,1)-sum(allFeatures,1))) < 1e-7, ...
        'Spatial block aggregation did not preserve global totals.');

    % Cluster bootstrap: sample nBlocks WHOLE geographic blocks WITH
    % replacement, using identical resample multiplicities for every
    % outcome. The number of pixels in a replicate may vary.
    rng(randomSeed+blockDeg,'twister');
    bootstrapSums = nan(nBootstrap,2*nFeature);
    for b = 1:nBootstrap
        selectedBlockID = randi(nBlocks,nBlocks,1);
        blockMultiplicity = accumarray( ...
            selectedBlockID,1,[nBlocks,1],@sum,0);
        bootstrapSums(b,:) = double(blockMultiplicity') * blockFeatures;
    end

    bootstrapResults(iScale).blockDegrees = blockDeg;
    bootstrapResults(iScale).nBlocks = nBlocks;
    bootstrapResults(iScale).blockSizes = blockSizes;
    bootstrapResults(iScale).blockID = blockID;

    pixelEffectiveCount = bootstrapSums(:,1);
    bootstrapResults(iScale).samplePixelCountMin = min(pixelEffectiveCount);
    bootstrapResults(iScale).samplePixelCountMax = max(pixelEffectiveCount);
    bootstrapResults(iScale).samplePixelCountMean = mean(pixelEffectiveCount);

    blockRows(end+1,:) = { ...
        blockDeg,nBlocks,min(blockSizes),median(blockSizes), ...
        max(blockSizes),mean(blockSizes), ...
        min(pixelEffectiveCount),max(pixelEffectiveCount)}; %#ok<SAGROW>

    for iWeight = 1:nWeighting
        modeName = weightingNames{iWeight};
        colStart = (iWeight-1)*nFeature+1;
        colEnd = iWeight*nFeature;
        sample = extractBootstrapStatistics( ...
            bootstrapSums(:,colStart:colEnd),rfxIndex,pairModelIndices);

        % Save all bootstrap replicate statistics for reproducibility.
        if iWeight == 1
            bootstrapResults(iScale).gridEqual = sample;
        else
            bootstrapResults(iScale).areaWeighted = sample;
        end
        point = pointResults(iWeight);

        for m = 1:nModel
            CI = percentileInterval(sample.modelMeans(:,m),confidenceLevel);
            meanRows(end+1,:) = { ...
                blockDeg,modeName,modelNamesExpected{m}, ...
                point.modelMeans(m),CI(1),CI(2),nBootstrap,nBlocks}; %#ok<SAGROW>
        end

        for p = 1:nPair
            % Standard pairwise 95% percentile interval.
            CI = percentileInterval( ...
                sample.pairedDifferences(:,p),confidenceLevel);
            excludesZero = CI(1)>0 || CI(2)<0;
            rfxLower = CI(2)<0;

            % Bonferroni simultaneous inference across all eight RFX
            % comparisons, within this block scale and weighting.
            % Same paired resamples: no additional bootstrap runs.
            CIB = percentileInterval( ...
                sample.pairedDifferences(:,p),bonferroniConfidenceLevel);
            bonfExcludesZero = CIB(1)>0 || CIB(2)<0;
            bonfRFXLower = CIB(2)<0;

            pairRows(end+1,:) = { ...
                blockDeg,modeName,pairNames{p}, ...
                point.pairedDifferences(p),CI(1),CI(2), ...
                excludesZero,rfxLower,nBootstrap,nBlocks, ...
                CIB(1),CIB(2),bonfExcludesZero,bonfRFXLower, ...
                bonferroniConfidenceLevel}; %#ok<SAGROW>
        end

        extendedCategoryNames = [categoryNames; ...
            {'Structural gain (Linear + Nonlinear)'}];
        extendedCategoryPoints = [point.categoryPct,point.structuralPct];
        extendedCategorySamples = [sample.categoryPct,sample.structuralPct];
        for c = 1:numel(extendedCategoryNames)
            CI = percentileInterval(extendedCategorySamples(:,c),confidenceLevel);
            categoryRows(end+1,:) = { ...
                blockDeg,modeName,extendedCategoryNames{c}, ...
                extendedCategoryPoints(c),CI(1),CI(2), ...
                nBootstrap,nBlocks}; %#ok<SAGROW>
        end

        for c = 1:nComponent
            CI = percentileInterval(sample.componentPct(:,c),confidenceLevel);
            componentRows(end+1,:) = { ...
                blockDeg,modeName,componentNames{c}, ...
                point.componentPct(c),CI(1),CI(2), ...
                nBootstrap,nBlocks}; %#ok<SAGROW>
        end

        for c = 1:2
            CI = percentileInterval(sample.contrasts(:,c),confidenceLevel);
            excludesZero = CI(1)>0 || CI(2)<0;
            contrastRows(end+1,:) = { ...
                blockDeg,modeName,contrastNames{c}, ...
                point.contrasts(c),CI(1),CI(2), ...
                excludesZero,nBootstrap,nBlocks}; %#ok<SAGROW>
        end

        fprintf('  %-12s | RFX-AR = %+.6f [%.6f, %.6f] | ', ...
            modeName,point.pairedDifferences(1), ...
            percentileIntervalValue(sample.pairedDifferences(:,1),confidenceLevel,1), ...
            percentileIntervalValue(sample.pairedDifferences(:,1),confidenceLevel,2));
        CIh = percentileInterval(sample.componentPct(:,4),confidenceLevel);
        fprintf('Hydro share = %.2f%% [%.2f, %.2f]\n', ...
            point.componentPct(4),CIh(1),CIh(2));
    end

    fprintf('  Geographic blocks: %d; median valid cells/block: %.1f\n\n', ...
        nBlocks,median(blockSizes));
    clear bootstrapSums blockFeatures;
end

%% ===================== 6. Build statistical tables ===================
modelCITable = cell2table(meanRows,'VariableNames', { ...
    'BlockDegrees','Weighting','Model','Estimate','CI_Lower', ...
    'CI_Upper','Replicates','NumBlocks'});
pairedCITable = cell2table(pairRows,'VariableNames', { ...
    'BlockDegrees','Weighting','Comparison','Estimate','CI_Lower', ...
    'CI_Upper','CI_ExcludesZero','CI_SupportsLowerRFX', ...
    'Replicates','NumBlocks', ...
    'BonferroniCI_Lower','BonferroniCI_Upper', ...
    'BonferroniCI_ExcludesZero','BonferroniCI_SupportsLowerRFX', ...
    'BonferroniIndividualConfidencePct'});

% The 3-degree, equal-grid table is the preselected main-paper comparison.
% All other block sizes and area weighting remain in the full CI table.
main3degPairTable = pairedCITable( ...
    pairedCITable.BlockDegrees==3 & ...
    strcmp(pairedCITable.Weighting,'GridEqual'),:);
assert(height(main3degPairTable)==nPair, ...
    'Main 3-degree comparison table must contain eight RFX contrasts.');
categoryCITable = cell2table(categoryRows,'VariableNames', { ...
    'BlockDegrees','Weighting','GainCategory','EstimatePct','CI_LowerPct', ...
    'CI_UpperPct','Replicates','NumBlocks'});
componentCITable = cell2table(componentRows,'VariableNames', { ...
    'BlockDegrees','Weighting','HydroClassComponent','EstimatePct', ...
    'CI_LowerPct','CI_UpperPct','Replicates','NumBlocks'});
contrastCITable = cell2table(contrastRows,'VariableNames', { ...
    'BlockDegrees','Weighting','Contrast','EstimatePP','CI_LowerPP', ...
    'CI_UpperPP','CI_ExcludesZero','Replicates','NumBlocks'});
blockTable = cell2table(blockRows,'VariableNames', { ...
    'BlockDegrees','NumBlocks','MinCellsPerBlock','MedianCellsPerBlock', ...
    'MaxCellsPerBlock','MeanCellsPerBlock', ...
    'MinResampledCells','MaxResampledCells'});

%% ========================== 7. Save outputs ==========================
% All reported point estimates retain their original statistical domains.
writetable(modelTable,fullfile(outDir,'s6_global_model_means.csv'));
writetable(categoryTable,fullfile(outDir,'s6_global_category_shares.csv'));
writetable(componentTable,fullfile(outDir,'s6_global_hydro_components.csv'));
writetable(modelCITable,fullfile(outDir,'s6_model_means_CI.csv'));
writetable(pairedCITable,fullfile(outDir,'s6_paired_model_differences_CI.csv'));
writetable(main3degPairTable,fullfile(outDir, ...
    's6_paired_model_differences_3deg.csv'));
writetable(categoryCITable,fullfile(outDir,'s6_gain_categories_CI.csv'));
writetable(componentCITable,fullfile(outDir,'s6_hydro_components_CI.csv'));
writetable(contrastCITable,fullfile(outDir,'s6_gain_contrasts_CI.csv'));
writetable(blockTable,fullfile(outDir,'s6_spatial_block_diagnostics.csv'));

excelFile = fullfile(outDir,'s6_spatial_block_bootstrap_tables.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(modelTable,excelFile,'Sheet','GlobalModels');
writetable(categoryTable,excelFile,'Sheet','GlobalCategories');
writetable(componentTable,excelFile,'Sheet','GlobalHydroComponents');
writetable(modelCITable,excelFile,'Sheet','ModelMeansCI');
writetable(pairedCITable,excelFile,'Sheet','PairedDiffCI');
writetable(main3degPairTable,excelFile,'Sheet','Paired3deg');
writetable(categoryCITable,excelFile,'Sheet','GainCategoriesCI');
writetable(componentCITable,excelFile,'Sheet','HydroComponentsCI');
writetable(contrastCITable,excelFile,'Sheet','GainContrastsCI');
writetable(blockTable,excelFile,'Sheet','SpatialBlocks');

%% ========================== 8. Optional figures ======================
if makeFigures
    try
        makeSummaryFigures(pairedCITable,categoryCITable, ...
            componentCITable,blockDegrees,outDir,settings);
        makeAllModelContrastFigure(main3degPairTable,outDir,settings);
    catch ME
        warning('Summary figures could not be produced: %s',ME.message);
    end
end

%% ============================ 9. Save MAT ============================
metadata.projectRoot = projectRoot;
metadata.stage1File = stage1File;
metadata.attributionFile = attributionFile;
metadata.blockDegrees = blockDegrees;
metadata.nBootstrap = nBootstrap;
metadata.confidenceLevel = confidenceLevel;
metadata.confidenceMethod = 'Two-sided percentile of cluster bootstrap';
metadata.multipleComparisonMethod = 'Bonferroni (eight RFX-vs-other contrasts)';
metadata.nPairedComparisons = nPair;
metadata.bonferroniIndividualConfidenceLevel = bonferroniConfidenceLevel;
metadata.pairModelIndices = pairModelIndices;
metadata.rfxIndex = rfxIndex;
metadata.randomSeed = randomSeed;
metadata.uncertaintyScope = 'Conditional spatial block resampling';
metadata.mainWeighting = 'GridEqual';
metadata.auxiliaryWeighting = 'AreaWeighted: cosd(latitude)';
metadata.testYears = [2022 2023 2024];

matFile = fullfile(outDir,'s6_spatial_block_bootstrap_results.mat');
save(matFile, ...
    'metadata','modelNamesExpected','categoryNames','componentNames', ...
    'pairNames','pairModelIndices','rfxIndex','contrastNames', ...
    'weightingNames','lat','lon', ...
    'validMask','validIndex','actualCategoryCounts', ...
    'pointResults','bootstrapResults','modelTable','categoryTable', ...
    'componentTable','modelCITable','pairedCITable','main3degPairTable', ...
    'categoryCITable','componentCITable','contrastCITable','blockTable', ...
    'maxReferenceError','maxARDifference','maxTotalGainDifference', ...
    'maxReconstructionError','settings','-v7.3');

fprintf('==================== GLOBAL POINT ESTIMATES ====================\n');
disp(modelTable);
disp(categoryTable);
disp(componentTable);
fprintf('================ MAIN 3-DEGREE RFX CONTRASTS ================\n');
disp(main3degPairTable(:,{ ...
    'Comparison','Estimate','CI_Lower','CI_Upper', ...
    'BonferroniCI_Lower','BonferroniCI_Upper', ...
    'BonferroniCI_SupportsLowerRFX'}));
fprintf('====================== BLOCK DIAGNOSTICS ======================\n');
disp(blockTable);
fprintf('Result MAT  : %s\n',matFile);
fprintf('Result XLSX : %s\n',excelFile);
fprintf('Elapsed time: %.1f seconds\n',toc(startTic));
fprintf('S6 COMPLETE. No forecasting model was retrained.\n');
fprintf('NOTE: CIs quantify spatial resampling under specified block sizes,\n');
fprintf('not new independent years, causal effects, or winner-selection bias.\n');
fprintf('===============================================================\n');
diary off;

%% =========================== Local functions ==========================
function checkMATVariables(fileName,requiredNames)
    if ~isfile(fileName)
        error('Required MAT file not found:\n%s',fileName);
    end
    existingNames = who('-file',fileName);
    missingNames = setdiff(requiredNames,existingNames);
    if ~isempty(missingNames)
        error('File %s is missing variable(s): %s', ...
            fileName,strjoin(missingNames,', '));
    end
end

function [blockID,nBlock] = buildSpatialBlocks(latCell,lonCell,blockDegrees)
% Fixed, globally aligned geographic rectangles. Edges are measured in
% degrees, NOT distance/area; all occupied blocks are sampled uniformly.
    latCell = double(latCell(:));
    lonCell = double(lonCell(:));
    assert(all(isfinite(latCell) & isfinite(lonCell)), ...
        'Valid-cell geographic coordinates must be finite.');
    assert(all(latCell>=-90 & latCell<=90), ...
        'Latitude values must lie between -90 and 90.');

    latBin = floor((latCell+90)./blockDegrees);
    lonWrapped = mod(lonCell+180,360);
    lonBin = floor(lonWrapped./blockDegrees);

    % unique(...,'rows') provides stable membership for all cells in a block.
    [~,~,blockID] = unique([latBin,lonBin],'rows');
    nBlock = max(blockID);
    assert(nBlock>1, 'Too few occupied blocks for cluster bootstrap.');
end

function result = extractStatistics(totals,rfxIndex,pairModelIndices)
% Reconstruct all metrics as ratios of SUFFICIENT-STATISTIC sums.
% For hydro components the denominator is the summed TOTAL gain in the
% selected Hydroclimatic-gain category, not the number of its grid cells.
    totals = double(totals(:)');
    denominator = totals(1);
    assert(isfinite(denominator) && denominator>0, ...
        'Cannot normalize a bootstrap sample with zero total weight.');
    result.modelMeans = totals(2:10)./denominator;
    result.categoryPct = 100.*totals(11:16)./denominator;
    result.structuralPct = sum(result.categoryPct(2:3));

    g = totals(17:20);
    assert(isfinite(g(4)) && g(4)>1e-12, ...
        'Hydro-class gain denominator is nonpositive in a replicate.');
    result.componentPct = 100.*[g(1),g(2),g(1)+g(2),g(3)]./g(4);
    result.pairedDifferences = ...
        result.modelMeans(rfxIndex)-result.modelMeans(pairModelIndices);
    result.contrasts = [ ...
        result.structuralPct-result.categoryPct(4), ...
        result.componentPct(3)-result.componentPct(4)];
end

function result = extractBootstrapStatistics(totals,rfxIndex,pairModelIndices)
% Vectorized equivalent of extractStatistics for nBootstrap rows.
    totals = double(totals);
    weight = totals(:,1);
    assert(all(isfinite(weight) & weight>0), ...
        'A resample has zero/nonfinite total weight.');
    result.modelMeans = bsxfun(@rdivide,totals(:,2:10),weight);
    result.categoryPct = 100.*bsxfun(@rdivide,totals(:,11:16),weight);
    result.structuralPct = sum(result.categoryPct(:,2:3),2);

    g = totals(:,17:20);
    assert(all(isfinite(g(:,4)) & g(:,4)>1e-12), ...
        'A bootstrap draw has nonpositive Hydro-class summed total gain.');
    result.componentPct = 100.*bsxfun(@rdivide, ...
        [g(:,1),g(:,2),g(:,1)+g(:,2),g(:,3)],g(:,4));
    result.pairedDifferences = bsxfun(@minus, ...
        result.modelMeans(:,rfxIndex), ...
        result.modelMeans(:,pairModelIndices));
    result.contrasts = [ ...
        result.structuralPct-result.categoryPct(:,4), ...
        result.componentPct(:,3)-result.componentPct(:,4)];
end

function CI = percentileInterval(values,confidenceLevel)
% A two-sided nonparametric PERCENTILE interval. No assumption that
% grid cells are independent; the resampling unit is the geographic block.
    values = double(values(:));
    assert(~isempty(values) && all(isfinite(values)), ...
        'Bootstrap confidence interval received nonfinite values.');
    alpha = (100-confidenceLevel)./2;
    CI = prctile(values,[alpha,100-alpha]);
end

function v = percentileIntervalValue(values,confidenceLevel,endpoint)
    CI = percentileInterval(values,confidenceLevel);
    v = CI(endpoint);
end

function makeSummaryFigures(pairedCI,categoryCI,componentCI, ...
    blockDegrees,outDir,settings)
% The diagrams use GRID-EQUAL results (main paper definition); area-weighted
% estimates remain available in the tables and MAT file.

    fig = figure('Color','w','Visible','off', ...
        'Position',settings.figPosition);
    pairLabels = {'RFX minus AR','RFX minus RF'};
    pairColors = [0.10 0.35 0.75; 0.75 0.30 0.15];

    for iPair=1:2
        subplot(1,2,iPair);
        rows = strcmp(pairedCI.Weighting,'GridEqual') & ...
            strcmp(pairedCI.Comparison,pairLabels{iPair});
        T = sortrows(pairedCI(rows,:),'BlockDegrees');
        x = T.BlockDegrees;
        y = T.Estimate;
        lower = y-T.CI_Lower;
        upper = T.CI_Upper-y;
        errorbar(x,y,lower,upper,'-o','LineWidth',1.4, ...
            'MarkerSize',6,'Color',pairColors(iPair,:));
        hold on;
        yline(0,'k--','LineWidth',0.8);
        xlim([min(blockDegrees)-0.4,max(blockDegrees)+0.4]);
        xticks(blockDegrees);
        xlabel('Spatial block size (degrees)');
        ylabel('Mean NRMSE difference');
        title(pairLabels{iPair},'Interpreter','none');
        set(gca,'FontName',settings.fontName, ...
            'FontSize',settings.tickFontSize);
        grid on;
    end
    exportgraphics(fig,fullfile(outDir, ...
        'Fig1_Paired_Model_Differences_CI.png'), ...
        'Resolution',settings.outputResolution);
    close(fig);

    fig = figure('Color','w','Visible','off', ...
        'Position',settings.figPosition);
    subplot(1,2,1);
    categoryLabels = { ...
        'Structural gain (Linear + Nonlinear)','Hydroclimatic gain'};
    categoryColors = [0.10 0.35 0.75; 0.10 0.60 0.30];
    for c = 1:2
        rows = strcmp(categoryCI.Weighting,'GridEqual') & ...
            strcmp(categoryCI.GainCategory,categoryLabels{c});
        T = sortrows(categoryCI(rows,:),'BlockDegrees');
        errorbar(T.BlockDegrees,T.EstimatePct, ...
            T.EstimatePct-T.CI_LowerPct, ...
            T.CI_UpperPct-T.EstimatePct, ...
            '-o','LineWidth',1.3,'MarkerSize',5, ...
            'Color',categoryColors(c,:));
        hold on;
    end
    xlabel('Spatial block size (degrees)');
    ylabel('Share of valid grid cells (%)');
    title('Global gain-category proportions');
    xticks(blockDegrees);
    xlim([min(blockDegrees)-0.4,max(blockDegrees)+0.4]);
    legend({'Linear + Nonlinear','Hydroclimatic'},'Location','best', ...
        'FontSize',settings.legendFontSize);
    set(gca,'FontName',settings.fontName, ...
        'FontSize',settings.tickFontSize);
    grid on;

    subplot(1,2,2);
    componentLabels = { ...
        'History (Linear + Nonlinear)','Hydroclimatic component'};
    for c=1:2
        rows = strcmp(componentCI.Weighting,'GridEqual') & ...
            strcmp(componentCI.HydroClassComponent,componentLabels{c});
        T = sortrows(componentCI(rows,:),'BlockDegrees');
        errorbar(T.BlockDegrees,T.EstimatePct, ...
            T.EstimatePct-T.CI_LowerPct, ...
            T.CI_UpperPct-T.EstimatePct, ...
            '-o','LineWidth',1.3,'MarkerSize',5, ...
            'Color',categoryColors(c,:));
        hold on;
    end
    xlabel('Spatial block size (degrees)');
    ylabel('Share of summed gain (%)');
    title('Within Hydroclimatic-gain category');
    xticks(blockDegrees);
    xlim([min(blockDegrees)-0.4,max(blockDegrees)+0.4]);
    legend({'History','Hydroclimatic'},'Location','best', ...
        'FontSize',settings.legendFontSize);
    set(gca,'FontName',settings.fontName, ...
        'FontSize',settings.tickFontSize);
    grid on;
    exportgraphics(fig,fullfile(outDir, ...
        'Fig2_Attribution_Uncertainty_CI.png'), ...
        'Resolution',settings.outputResolution);
    close(fig);
end

function makeAllModelContrastFigure(T,outDir,settings)
% Main 3-degree, grid-equal eight-model comparison with simultaneous CIs.
% Conventional 95% intervals are reported in tables; broader Bonferroni
% intervals are plotted here to prevent selective significance claims.
    assert(height(T)==8 && all(T.BlockDegrees==3) && ...
        all(strcmp(T.Weighting,'GridEqual')), ...
        'The all-model figure requires eight 3-degree grid-equal rows.');
    otherNames = strrep(T.Comparison,'RFX minus ','');
    f = figure('Color','w','Visible','off', ...
        'Position',[100 100 1100 700]);
    ax = axes('Parent',f);
    hold(ax,'on');
    n = height(T);
    y = (1:n)';
    negativeColor = [0.10 0.40 0.68];
    overlapColor = [0.66 0.34 0.23];
    for k=1:n
        if T.BonferroniCI_SupportsLowerRFX(k)
            lineColor = negativeColor;
        else
            lineColor = overlapColor;
        end
        plot(ax,[T.BonferroniCI_Lower(k),T.BonferroniCI_Upper(k)], ...
            [y(k),y(k)],'-','Color',lineColor,'LineWidth',2);
        plot(ax,[T.BonferroniCI_Lower(k),T.BonferroniCI_Lower(k)], ...
            [y(k)-0.12,y(k)+0.12],'-','Color',lineColor);
        plot(ax,[T.BonferroniCI_Upper(k),T.BonferroniCI_Upper(k)], ...
            [y(k)-0.12,y(k)+0.12],'-','Color',lineColor);
        plot(ax,T.Estimate(k),y(k),'o','MarkerSize',6, ...
            'MarkerFaceColor',lineColor,'MarkerEdgeColor','w');
    end
    xline(ax,0,'k--','LineWidth',0.9);
    set(ax,'YTick',1:n,'YTickLabel',otherNames, ...
        'YDir','reverse','FontName',settings.fontName, ...
        'FontSize',settings.tickFontSize);
    ylim(ax,[0.4,n+0.6]);
    xlabel(ax,'Global mean NRMSE difference (RFX minus comparison model)', ...
        'FontName',settings.fontName, ...
        'FontSize',settings.axisLabelFontSize);
    ylabel(ax,'Comparison model','FontName',settings.fontName, ...
        'FontSize',settings.axisLabelFontSize);
    title(ax,'3-degree spatial blocks: Bonferroni-adjusted CIs (8 contrasts)', ...
        'FontName',settings.fontName,'FontWeight','normal');
    grid(ax,'on');
    exportgraphics(f,fullfile(outDir, ...
        'Fig3_RFX_All_Model_Contrasts_3deg.png'), ...
        'Resolution',settings.outputResolution);
    close(f);
end
