%% STEP 4 - Upscale GLDAS-DA groundwater storage from 0.25 to 1 degree
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Aggregate GLDAS-DA V2.2 monthly-end GWS to a 1-degree grid and calculate
%   monthly groundwater storage changes from adjacent monthly-end states.
%
% Input
%   GLDAS_DA_CLSM025_GWS_monthend_20031231_20241231.mat (Step 3).
%   GWS_mm: [lat, lon, state], 2003-12-31 through 2024-12-31 (253 dates).
%   Each GWS_tavg record is the daily mean on the final day of its month.
%
% Spatial aggregation
%   Each 1-degree cell contains 4 x 4 native 0.25-degree cells.
%   Valid native values are weighted by spherical cell area:
%       X_1deg = sum(A_j * X_j) / sum(A_j), over valid native cells.
%   A_j is proportional to sind(lat_j + 0.125) - sind(lat_j - 0.125).
%   Longitude width and Earth-radius factors cancel in the weighted mean.
%   Output is NaN when all native cells in a block are invalid.
%
% Monthly change
%   Calculate native-cell changes before area-weighted aggregation:
%       dGWS025(t) = GWS025(t) - GWS025(t-1).
%   Only native cells valid in BOTH months contribute to a monthly change.
%   This also defines how monthly soil-water changes are aggregated for Noah.
%   GWS1_mm is independently aggregated for each monthly-end state. When
%   source-cell coverage varies between months, dGWS_mm need not equal
%   diff(GWS1_mm, 1, 3); the difference is reported as a diagnostic.
%
% Output
%   GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat
%   GWS1_mm: [lat, lon, state], 253 dates; unit: mm.
%   dGWS_mm: [lat, lon, time], 2004-01 to 2024-12; unit: mm/month.
%   lat, lon, Lat2D, Lon2D: ascending 1-degree grid coordinates.
%   ym, yearVec, monthVec: 252 monthly change dates.
%   ymState, yyyymmddState: 253 monthly-end state dates.
%   ValidAreaFraction.GWS1_mm: source area valid in each state (0 to 1).
%   ValidAreaFraction.dGWS_mm: source area valid in both months (0 to 1).
%   Meta: product information and calculation definitions. MAT format: -v7.3.
%
% Reference
%   NASA GES DISC, README_GLDAS2.pdf, Sections 3.1 and 3.5.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
clearvars;
clc;

projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised\s1_master_database';
inFile = fullfile(projectRoot, '01_GLDAS_extraction', ...
    'GLDAS_DA_CLSM025_GWS_monthend_20031231_20241231.mat');
outFolder = fullfile(projectRoot, '02_spatial_temporal_processing');
outFileName = 'GLDAS_DA_CLSM10_dGWS_monthend_200401_202412.mat';

% Existing output file: "overwrite" or "stop".
outputMode = "overwrite";

% Number of monthly changes processed per block.
chunkMonths = 12;

makeDiagnosticPlot = true;
plotYM = 200401;

%% =========================================================================
% 3. MAIN EXECUTION
% =========================================================================

%% 3.1 Check paths and parameters
if ~isfile(inFile)
    error('Input MAT-file not found: %s', inFile);
end
if ~ismember(outputMode, ["overwrite", "stop"])
    error('outputMode must be "overwrite" or "stop".');
end
if ~isscalar(chunkMonths) || ~isfinite(chunkMonths) || ...
        chunkMonths < 1 || chunkMonths ~= floor(chunkMonths)
    error('chunkMonths must be a positive integer.');
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

%% 3.2 Read coordinates and monthly-end dates
requiredVars = {'GWS_mm', 'lat', 'lon', 'ymState', 'yyyymmddState'};
fileVars = who('-file', inFile);
for k = 1:numel(requiredVars)
    if ~any(strcmp(fileVars, requiredVars{k}))
        error('Input MAT-file is missing %s.', requiredVars{k});
    end
end

S = matfile(inFile);
Coord = load(inFile, 'lat', 'lon', 'ymState', 'yyyymmddState');
lat025 = double(Coord.lat(:));
lon025 = double(Coord.lon(:));
ymState = double(Coord.ymState(:));
yyyymmddState = double(Coord.yyyymmddState(:));
clear Coord

validateNativeGrid(lat025, lon025);
nLat025 = numel(lat025);
nLon025 = numel(lon025);
nLat1 = nLat025 / 4;
nLon1 = nLon025 / 4;

expectedYM = makeYMList(2003, 12, 2024, 12);
if ~isequal(ymState, expectedYM) || numel(yyyymmddState) ~= numel(ymState)
    error('Input must contain 253 consecutive states from 2003-12 to 2024-12.');
end
expectedDates = floor(ymState / 100) * 10000 + ...
    mod(ymState, 100) * 100 + ...
    eomday(floor(ymState / 100), mod(ymState, 100));
if ~isequal(yyyymmddState, expectedDates)
    error('State dates do not match the final calendar days of the months.');
end

expectedNativeSize = [nLat025, nLon025, numel(ymState)];
actualSize = size(S, 'GWS_mm');
if ~isequal(actualSize, expectedNativeSize)
    error('GWS_mm has size %s; expected %s.', ...
        mat2str(actualSize), mat2str(expectedNativeSize));
end

ym = ymState(2:end);
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);
nTime = numel(ym);
nState = numel(ymState);

fprintf('Input file: %s\n', inFile);
fprintf('Native grid: %d latitude x %d longitude x %d state dates\n', ...
    nLat025, nLon025, nState);
fprintf('Target period: %d to %d (%d months)\n', ym(1), ym(end), nTime);

%% 3.3 Define 1-degree coordinates and spherical area weights
lat = mean(reshape(lat025, 4, nLat1), 1)';
lon = mean(reshape(lon025, 4, nLon1), 1)';
[Lon2D, Lat2D] = meshgrid(lon, lat);

rowAreaWeight = sind(lat025 + 0.125) - sind(lat025 - 0.125);
if any(~isfinite(rowAreaWeight)) || any(rowAreaWeight <= 0)
    error('Invalid spherical area weights for the source latitude grid.');
end

fprintf('Target grid: %d latitude x %d longitude (1 degree)\n', nLat1, nLon1);
fprintf('Using spherical cell-area weights with available-area normalization.\n');

%% 3.4 Allocate output arrays
GWS1_mm = NaN(nLat1, nLon1, nState);
dGWS_mm = NaN(nLat1, nLon1, nTime);
ValidAreaFraction = struct();
ValidAreaFraction.GWS1_mm = NaN(nLat1, nLon1, nState, 'single');
ValidAreaFraction.dGWS_mm = NaN(nLat1, nLon1, nTime, 'single');

%% 3.5 Aggregate monthly-end GWS and monthly changes
runTimer = tic;
for t1 = 1:chunkMonths:nTime
    t2 = min(t1 + chunkMonths - 1, nTime);

    % One reference state is included for each block of monthly changes.
    GWS025 = double(S.GWS_mm(:, :, t1:(t2+1)));
    [stateBlock, stateCoverage] = areaWeighted025to1(GWS025, rowAreaWeight);

    if t1 == 1
        GWS1_mm(:, :, 1:(t2+1)) = stateBlock;
        ValidAreaFraction.GWS1_mm(:, :, 1:(t2+1)) = stateCoverage;
    else
        % The first state in this block was stored by the preceding block.
        GWS1_mm(:, :, (t1+1):(t2+1)) = stateBlock(:, :, 2:end);
        ValidAreaFraction.GWS1_mm(:, :, (t1+1):(t2+1)) = ...
            stateCoverage(:, :, 2:end);
    end

    nativeChange = diff(GWS025, 1, 3);
    [changeBlock, changeCoverage] = areaWeighted025to1( ...
        nativeChange, rowAreaWeight);
    dGWS_mm(:, :, t1:t2) = changeBlock;
    ValidAreaFraction.dGWS_mm(:, :, t1:t2) = changeCoverage;

    fprintf('Processed %3d/%3d monthly changes\n', t2, nTime);
end
fprintf('Aggregation completed in %.2f minutes.\n', toc(runTimer) / 60);

%% 3.6 Check output arrays and temporal support
checkOutputSize(GWS1_mm, [nLat1, nLon1, nState], 'GWS1_mm');
checkOutputSize(dGWS_mm, [nLat1, nLon1, nTime], 'dGWS_mm');

fprintf('\nOutput variable summary:\n');
printVarSummary(GWS1_mm, ValidAreaFraction.GWS1_mm, 'GWS1_mm');
printVarSummary(dGWS_mm, ValidAreaFraction.dGWS_mm, 'dGWS_mm');

% Differences between pairwise-valid changes and differences of monthly
% means indicate changes in native-cell availability between state dates.
stateMeanChange = diff(GWS1_mm, 1, 3);
commonFinite = isfinite(dGWS_mm) & isfinite(stateMeanChange);
if any(commonFinite(:))
    absDifference = abs(dGWS_mm(commonFinite) - stateMeanChange(commonFinite));
    supportDifferenceCount = nnz(absDifference > 1e-8);
    fprintf(['Native-paired vs. state-mean dGWS: %.6g mm mean absolute ' ...
        'difference; %d / %d cases exceed 1e-8 mm.\n'], ...
        mean(absDifference), supportDifferenceCount, numel(absDifference));
else
    fprintf('No paired valid changes available for the temporal comparison.\n');
end
clear stateMeanChange commonFinite absDifference

%% 3.7 Metadata and output
Meta = struct();
Meta.sourceProduct = 'GLDAS_CLSM025_DA1_D_2.2';
Meta.sourceDOI = '10.5067/TXBMLX370XX8';
Meta.periodState = 'Monthly-end daily-mean GWS, 2003-12-31 to 2024-12-31';
Meta.periodTarget = 'Monthly dGWS, 2004-01 to 2024-12';
Meta.arrayConvention = '[lat, lon, time], ascending latitude and longitude';
Meta.spatialUpscaling = ['4 x 4 native cells per 1-degree cell; ' ...
    'spherical cell-area-weighted mean over finite values'];
Meta.areaWeight = 'sind(lat + 0.125) - sind(lat - 0.125)';
Meta.missingDataRule = 'NaN when no valid source cells; otherwise normalize by valid area';
Meta.ValidAreaFraction = ['GWS1_mm: valid area at each state; dGWS_mm: ' ...
    'area with valid GWS on both adjacent month-end dates'];
Meta.GWS1_mm = 'Area-weighted daily-mean GWS on each selected month-end date; mm';
Meta.dGWS_mm = ['GWS025(t) - GWS025(t-1), then area-weighted over ' ...
    'native cells valid in both months; mm/month'];
Meta.temporalSupport = ['Independent state averages and paired native-cell ' ...
    'changes may differ if the valid native-cell mask changes.'];
Meta.targetExample1 = '2004-01 dGWS = GWS(2004-01-31) - GWS(2003-12-31)';
Meta.targetExample2 = '2024-12 dGWS = GWS(2024-12-31) - GWS(2024-11-30)';

fprintf('\nSaving MAT-file: %s\n', outFile);
save(outFile, ...
    'GWS1_mm', 'dGWS_mm', ...
    'lat', 'lon', 'Lat2D', 'Lon2D', ...
    'ym', 'yearVec', 'monthVec', 'ymState', 'yyyymmddState', ...
    'ValidAreaFraction', 'Meta', '-v7.3');
fprintf('Saved successfully.\n');

%% 3.8 Diagnostic map
if makeDiagnosticPlot
    plotIndex = find(ym == plotYM, 1);
    if isempty(plotIndex)
        warning('Diagnostic month %06d is outside the target period.', plotYM);
    else
        gwsMap = dGWS_mm(:, :, plotIndex);
        figure('Name', 'GLDAS-DA 1-degree monthly GWS change');
        h = imagesc(lon, lat, gwsMap);
        set(h, 'AlphaData', double(isfinite(gwsMap)));
        set(gca, 'YDir', 'normal');
        axis equal tight;
        box on;
        grid on;
        colormap(jet);
        cb = colorbar;
        ylabel(cb, 'dGWS (mm/month)');
        title(sprintf('GLDAS-DA CLSM monthly GWS change, %04d-%02d, 1 degree', ...
            floor(plotYM / 100), mod(plotYM, 100)));
        xlabel('Longitude');
        ylabel('Latitude');
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
    end
end

%% =========================================================================
% 4. LOCAL FUNCTIONS
% =========================================================================

function ym = makeYMList(y1, m1, y2, m2)
    d1 = datetime(y1, m1, 1);
    d2 = datetime(y2, m2, 1);
    if d2 < d1
        error('End month must not precede start month.');
    end
    d = d1:calmonths(1):d2;
    ym = double(year(d(:)) * 100 + month(d(:)));
end

function validateNativeGrid(lat025, lon025)
    if numel(lat025) ~= 600 || numel(lon025) ~= 1440 || ...
            any(~isfinite(lat025)) || any(~isfinite(lon025))
        error('Input must contain the full GLDAS 600 x 1440 native grid.');
    end
    tol = 1e-5;
    if any(abs(diff(lat025) - 0.25) > tol) || ...
            any(abs(diff(lon025) - 0.25) > tol) || ...
            abs(lat025(1) + 59.875) > tol || ...
            abs(lat025(end) - 89.875) > tol || ...
            abs(lon025(1) + 179.875) > tol || ...
            abs(lon025(end) - 179.875) > tol
        error('Input coordinates are not aligned to the GLDAS 0.25-degree grid.');
    end
end

function [Y, C] = areaWeighted025to1(X, wLat)
    % Input: [lat025, lon025, time]. Output: [lat1, lon1, time].
    [nLat025, nLon025, nTime] = size(X);
    if mod(nLat025, 4) ~= 0 || mod(nLon025, 4) ~= 0 || ...
            numel(wLat) ~= nLat025
        error('Incompatible source dimensions or area weights.');
    end
    nLat1 = nLat025 / 4;
    nLon1 = nLon025 / 4;

    % Reshape into [4 native latitudes, 1-degree latitudes,
    % 4 native longitudes, 1-degree longitudes, months].
    X = reshape(X, 4, nLat1, 4, nLon1, nTime);
    w = reshape(wLat, 4, nLat1, 1, 1, 1);
    valid = isfinite(X);
    X(~valid) = 0;

    weightedSum = sum(sum(X .* w, 1), 3);
    validWeight = sum(sum(double(valid) .* w, 1), 3);
    numerator = reshape(weightedSum, nLat1, nLon1, nTime);
    denominator = reshape(validWeight, nLat1, nLon1, nTime);

    Y = numerator ./ denominator;
    Y(denominator == 0) = NaN;

    totalWeight = 4 * sum(reshape(wLat, 4, nLat1), 1);
    C = single(denominator ./ reshape(totalWeight, nLat1, 1, 1));
    C = max(0, min(1, C));
end


function checkOutputSize(X, expectedSize, varName)
    if ~isequal(size(X), expectedSize) || ~isa(X, 'double')
        error('Unexpected size or data type for %s.', varName);
    end
end

function printVarSummary(X, C, varName)
    nFinite = nnz(isfinite(X));
    nPartial = nnz(C(:) > 0 & C(:) < 1 - 1e-6);
    fprintf(['%-10s | valid: %7.2f%% | partial-area cells: %8d | ' ...
        'min: %10.4f | mean: %10.4f | max: %10.4f\n'], ...
        varName, 100 * nFinite / numel(X), nPartial, ...
        min(X(:), [], 'omitnan'), mean(X(:), 'omitnan'), ...
        max(X(:), [], 'omitnan'));
end
