%% STEP 5 - Build the common GLDAS-DA and GLDAS-Noah valid-grid database
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Identify 1-degree grid cells with complete monthly groundwater-storage
%   change and hydroclimatic predictor records for 2004-01 to 2024-12.
%
% Input
%   GLDAS_NOAH10_5vars_200401_202412.mat (Step 2)
%       Prec_mm, Evap_mm, Runoff_mm, dSM2m_mm, SurfT_C
%   GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat (Step 4)
%       dGWS_mm
%   Each input field has dimensions [lat, lon, time] = [150, 360, 252].
%
% Mask rule
%   A grid cell is retained only when all six fields are finite in every
%   month. No minimum native-cell area-coverage threshold is imposed.
%   Noah and DA may have different native-cell coverage within a 1-degree
%   cell. Their coverage fractions are retained for data auditing.
%
% Output
%   GLDAS_common_valid_mask_1deg_200401_202412.mat
%   CommonMask:         [150, 360], 1 for valid cells and NaN otherwise
%   CommonValid:        [150, 360], logical common validity mask
%   validLinearIndex:   linear indices into [lat, lon] grids
%   validLat, validLon: coordinates of retained grid cells
%   dGWS_valid, Prec_valid, Evap_valid, Runoff_valid,
%   dSM2m_valid, SurfT_valid: [validGridCount, 252]
%   CoverageMin, CoverageMean: native-cell area-coverage fractions for the
%       six variables, summarized across 252 months at retained grid cells
%   MaskAudit: counts of cells satisfying each completeness requirement
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
clearvars;
clc;

projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised\s1_master_database';
inputFolder = fullfile(projectRoot, '02_spatial_temporal_processing');
noahFile = fullfile(inputFolder, 'GLDAS_NOAH10_5vars_200401_202412.mat');
daFile = fullfile(inputFolder, 'GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat');
outFolder = fullfile(projectRoot, '03_master_database');
outFileName = 'GLDAS_common_valid_mask_1deg_200401_202412.mat';

% Existing output file: "overwrite" or "stop".
outputMode = "overwrite";

% Manuscript grid count is checked for reference, not enforced.
reportedGridCount = 17284;

coordTol = 1e-8;
makeDiagnosticPlot = true;

%% =========================================================================
% 3. MAIN EXECUTION
% =========================================================================

%% 3.1 Check input files and output settings
if ~isfile(noahFile)
    error('Noah MAT-file not found: %s', noahFile);
end
if ~isfile(daFile)
    error('GLDAS-DA MAT-file not found: %s', daFile);
end
if ~ismember(outputMode, ["overwrite", "stop"])
    error('outputMode must be "overwrite" or "stop".');
end
if ~isfolder(outFolder)
    [ok, msg] = mkdir(outFolder);
    if ~ok
        error('Cannot create output directory: %s\n%s', outFolder, msg);
    end
end
outFile = fullfile(outFolder, outFileName);
if isfile(outFile) && outputMode == "stop"
    error('Output MAT-file already exists: %s', outFile);
end

noahVars = who('-file', noahFile);
daVars = who('-file', daFile);
requiredNoahVars = {'Prec_mm', 'Evap_mm', 'Runoff_mm', ...
    'dSM2m_mm', 'SurfT_C', 'lat', 'lon', 'ym', 'ValidAreaFraction'};
requiredDaVars = {'dGWS_mm', 'lat', 'lon', 'ym', 'ValidAreaFraction'};
checkVarsExist(noahVars, requiredNoahVars, 'GLDAS-Noah');
checkVarsExist(daVars, requiredDaVars, 'GLDAS-DA');

%% 3.2 Check coordinates and monthly sequence
Nmeta = load(noahFile, 'lat', 'lon', 'ym');
Dmeta = load(daFile, 'lat', 'lon', 'ym');

latNoah = double(Nmeta.lat(:));
lonNoah = double(Nmeta.lon(:));
ymNoah = double(Nmeta.ym(:));
latDA = double(Dmeta.lat(:));
lonDA = double(Dmeta.lon(:));
ymDA = double(Dmeta.ym(:));
clear Nmeta Dmeta

validateTargetGrid(latNoah, lonNoah, coordTol);
validateTargetGrid(latDA, lonDA, coordTol);
if any(abs(latNoah - latDA) > coordTol) || ...
        any(abs(lonNoah - lonDA) > coordTol)
    error('The GLDAS-Noah and GLDAS-DA grids are not aligned.');
end

expectedYM = makeYMList(2004, 1, 2024, 12);
if ~isequal(ymNoah, expectedYM) || ~isequal(ymDA, expectedYM)
    error('Both inputs must contain every month from 2004-01 to 2024-12.');
end

lat = latNoah;
lon = lonNoah;
ym = expectedYM;
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);
nLat = numel(lat);
nLon = numel(lon);
nTime = numel(ym);
[Lon2D, Lat2D] = meshgrid(lon, lat);
expectedSize = [nLat, nLon, nTime];

fprintf('Grid: %d latitudes x %d longitudes\n', nLat, nLon);
fprintf('Monthly period: %06d to %06d (%d months)\n', ym(1), ym(end), nTime);

%% 3.3 Load and validate six monthly fields
fprintf('\nLoading GLDAS-DA and GLDAS-Noah fields...\n');
D = load(daFile, 'dGWS_mm');
dGWS_mm = double(D.dGWS_mm);
clear D

N = load(noahFile, 'Prec_mm', 'Evap_mm', 'Runoff_mm', ...
    'dSM2m_mm', 'SurfT_C');
Prec_mm = double(N.Prec_mm);
Evap_mm = double(N.Evap_mm);
Runoff_mm = double(N.Runoff_mm);
dSM2m_mm = double(N.dSM2m_mm);
SurfT_C = double(N.SurfT_C);
clear N

fieldNames = {'dGWS_mm', 'Prec_mm', 'Evap_mm', ...
    'Runoff_mm', 'dSM2m_mm', 'SurfT_C'};
fields = {dGWS_mm, Prec_mm, Evap_mm, Runoff_mm, dSM2m_mm, SurfT_C};

for k = 1:numel(fieldNames)
    checkArraySize(fields{k}, expectedSize, fieldNames{k});
end

%% 3.4 Build the common full-period validity mask
fprintf('Checking the six variables over all %d months...\n', nTime);
MaskAudit = struct();
CommonValid = true(nLat, nLon);
for k = 1:numel(fieldNames)
    validThis = all(isfinite(fields{k}), 3);
    MaskAudit.validByVariable.(fieldNames{k}) = nnz(validThis);
    CommonValid = CommonValid & validThis;
    fprintf('%-12s: %6d complete grid cells\n', ...
        fieldNames{k}, nnz(validThis));
end

CommonMask = NaN(nLat, nLon);
CommonMask(CommonValid) = 1;
validLinearIndex = find(CommonValid);
validGridCount = numel(validLinearIndex);
totalGridCount = nLat * nLon;
validGridFraction = validGridCount / totalGridCount;

if validGridCount == 0
    error('No grid cells satisfy the six-variable completeness mask.');
end
fprintf('\nTotal 1-degree cells: %d\n', totalGridCount);
fprintf('Common valid cells : %d (%.2f%%)\n', ...
    validGridCount, 100 * validGridFraction);
if validGridCount ~= reportedGridCount
    warning('SKLGB:GridCount', ...
        ['The retained count (%d) differs from the manuscript count (%d). ' ...
         'The calculated mask is retained.'], validGridCount, reportedGridCount);
end

MaskAudit.totalGridCount = totalGridCount;
MaskAudit.commonValidGridCount = validGridCount;
MaskAudit.periodMonths = nTime;
MaskAudit.reportedGridCount = reportedGridCount;

%% 3.5 Extract monthly time series for the retained grid cells
validLat = Lat2D(validLinearIndex);
validLon = Lon2D(validLinearIndex);

dGWS_valid = extractValidSeries(dGWS_mm, validLinearIndex);
Prec_valid = extractValidSeries(Prec_mm, validLinearIndex);
Evap_valid = extractValidSeries(Evap_mm, validLinearIndex);
Runoff_valid = extractValidSeries(Runoff_mm, validLinearIndex);
dSM2m_valid = extractValidSeries(dSM2m_mm, validLinearIndex);
SurfT_valid = extractValidSeries(SurfT_C, validLinearIndex);

validFields = {dGWS_valid, Prec_valid, Evap_valid, ...
    Runoff_valid, dSM2m_valid, SurfT_valid};
expectedValidSize = [validGridCount, nTime];
for k = 1:numel(fieldNames)
    if ~isequal(size(validFields{k}), expectedValidSize) || ...
            any(~isfinite(validFields{k}(:)))
        error('Invalid extracted time-series matrix: %s', fieldNames{k});
    end
end
fprintf('\nExtracted six time-series matrices: [%d, %d]\n', ...
    validGridCount, nTime);

%% 3.6 Summarize native-cell area coverage
% Coverage summaries are diagnostic only and do not change CommonValid.
Ncov = load(noahFile, 'ValidAreaFraction');
Dcov = load(daFile, 'ValidAreaFraction');
noahCoverage = Ncov.ValidAreaFraction;
daCoverage = Dcov.ValidAreaFraction;
clear Ncov Dcov

CoverageMin = struct();
CoverageMean = struct();
for k = 1:numel(fieldNames)
    varName = fieldNames{k};
    if k == 1
        if ~isfield(daCoverage, varName)
            error('DA ValidAreaFraction is missing %s.', varName);
        end
        C = daCoverage.(varName);
    else
        if ~isfield(noahCoverage, varName)
            error('Noah ValidAreaFraction is missing %s.', varName);
        end
        C = noahCoverage.(varName);
    end
    checkCoverageConsistency(C, fields{k}, expectedSize, varName);
    C2D = reshape(C, nLat * nLon, nTime);
    Cvalid = C2D(validLinearIndex, :);
    CoverageMin.(varName) = min(Cvalid, [], 2);
    CoverageMean.(varName) = mean(Cvalid, 2);
    fprintf('%-12s: mean native coverage %.2f%%; partial cells %.2f%%\n', ...
        varName, 100 * mean(Cvalid(:)), ...
        100 * nnz(Cvalid(:) < 1 - 1e-6) / numel(Cvalid));
    clear C C2D Cvalid
end
clear noahCoverage daCoverage

%% 3.7 Save the modeling database
Meta = struct();
Meta.sourceProducts = 'GLDAS_CLSM025_DA1_D_2.2 and GLDAS_NOAH025_M_2.1';
Meta.period = '2004-01 to 2024-12';
Meta.arrayConvention = 'Spatial [lat, lon]; time series [grid, month]';
Meta.maskDefinition = ['Retain a 1-degree grid cell when all six variables ' ...
    'are finite for all 252 months; no minimum source-area threshold'];
Meta.spatialUpscaling = ['Both products use spherical source-cell area ' ...
    'weights normalized by the area of valid 0.25-degree source cells'];
Meta.coverage = ['CoverageMin and CoverageMean summarize each variable''s ' ...
    'native-cell area fraction over the retained monthly records'];
Meta.target = 'dGWS_valid from GLDAS-DA; mm/month';
Meta.predictors = ['Prec_valid, Evap_valid, Runoff_valid, dSM2m_valid, ' ...
    'SurfT_valid from GLDAS-Noah'];
Meta.scriptName = 'step5_create_mask.m';

fprintf('\nSaving MAT-file: %s\n', outFile);
save(outFile, ...
    'CommonMask', 'CommonValid', ...
    'validGridCount', 'totalGridCount', 'validGridFraction', ...
    'validLinearIndex', 'validLat', 'validLon', ...
    'dGWS_valid', 'Prec_valid', 'Evap_valid', ...
    'Runoff_valid', 'dSM2m_valid', 'SurfT_valid', ...
    'lat', 'lon', 'Lat2D', 'Lon2D', ...
    'ym', 'yearVec', 'monthVec', ...
    'CoverageMin', 'CoverageMean', 'MaskAudit', 'Meta', '-v7.3');
fprintf('Saved successfully.\n');

%% 3.8 Diagnostic map
if makeDiagnosticPlot
    figure('Name', 'GLDAS common valid-grid mask');
    imagesc(lon, lat, double(CommonValid));
    set(gca, 'YDir', 'normal');
    axis equal tight;
    box on;
    grid on;
    colormap(gray(2));
    caxis([0 1]);
    cb = colorbar;
    ylabel(cb, 'Validity (0 = excluded, 1 = retained)');
    title(sprintf('Common valid grid cells, n = %d', validGridCount));
    xlabel('Longitude');
    ylabel('Latitude');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

function ym = makeYMList(y1, m1, y2, m2)
    d1 = datetime(y1, m1, 1);
    d2 = datetime(y2, m2, 1);
    d = d1:calmonths(1):d2;
    ym = double(year(d(:)) * 100 + month(d(:)));
end

function validateTargetGrid(lat, lon, tol)
    if numel(lat) ~= 150 || numel(lon) ~= 360 || ...
            any(~isfinite(lat)) || any(~isfinite(lon)) || ...
            any(abs(diff(lat) - 1) > tol) || ...
            any(abs(diff(lon) - 1) > tol) || ...
            abs(lat(1) + 59.5) > tol || ...
            abs(lat(end) - 89.5) > tol || ...
            abs(lon(1) + 179.5) > tol || ...
            abs(lon(end) - 179.5) > tol
        error('Expected the full aligned GLDAS 150 x 360, 1-degree grid.');
    end
end

function checkVarsExist(existingVars, requiredVars, fileLabel)
    for k = 1:numel(requiredVars)
        if ~any(strcmp(existingVars, requiredVars{k}))
            error('%s is missing required variable %s.', fileLabel, requiredVars{k});
        end
    end
end

function checkArraySize(X, expectedSize, varName)
    if ~isequal(size(X), expectedSize)
        error('%s has size %s; expected %s.', ...
            varName, mat2str(size(X)), mat2str(expectedSize));
    end
end

function Xvalid = extractValidSeries(X, idx)
    [nLat, nLon, nTime] = size(X);
    X2D = reshape(X, nLat * nLon, nTime);
    Xvalid = X2D(idx, :);
end

function checkCoverageConsistency(C, X, expectedSize, varName)
    if ~isequal(size(C), expectedSize) || ...
            any(~isfinite(C(:))) || ...
            any(C(:) < -1e-6) || any(C(:) > 1 + 1e-6)
        error('Invalid native-cell coverage array for %s.', varName);
    end
    if any(isfinite(X(:)) ~= (C(:) > 0))
        error('Coverage and finite-value masks disagree for %s.', varName);
    end
end
