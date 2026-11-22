%% STEP 3 - Extract GLDAS-DA CLSM monthly-end groundwater storage
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Extract GWS_tavg from the final calendar day of each month in the
%   GLDAS-2.2 GRACE-DA1 Catchment LSM daily 0.25-degree product.
%
% Input product
%   GLDAS_CLSM025_DA1_D_2.2 (Main Production)
%   Files: GLDAS_CLSM025_DA1_D.AYYYYMMDD.022.nc4
%   NASA product DOI: 10.5067/TXBMLX370XX8
%
% Period
%   2003-12-31 to 2024-12-31 (253 monthly-end dates).
%   December 2003 supplies the reference state for January 2004 dGWS.
%
% Variable
%   GWS_mm = GWS_tavg [mm]. GWS_tavg is the daily mean groundwater
%   storage for the selected day, not an instantaneous month-end value.
%   Monthly storage changes are calculated in the spatial processing step.
%
% Output
%   GLDAS_DA_CLSM025_GWS_monthend_20031231_20241231.mat
%   GWS_mm [lat, lon, time]; lat south-to-north, lon west-to-east.
%   lat, lon, Lat2D, Lon2D: native 0.25-degree coordinates.
%   dateState, ymState, yyyymmddState: state dates.
%   yearVecState, monthVecState, dayVecState: date components.
%   Meta: source and processing information. MAT format: -v7.3.
%
% Reference
%   NASA GES DISC, README_GLDAS2.pdf, Sections 3.1 and 3.5.
%
%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
clearvars;
clc;

% Directory with daily GLDAS-DA NetCDF files arranged in year folders.
folderPath = ...
    'D:\A1_科研数据库\01_原始数据\01_GLDAS\GLDAS-DA CLSM Daily 025D 20030201-to-20251231';

% Output directory.
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised\s1_master_database';
outFolder = fullfile(projectRoot, '01_GLDAS_extraction');
outFileName = 'GLDAS_DA_CLSM025_GWS_monthend_20031231_20241231.mat';

% Output file handling: "overwrite" replaces the file; "stop" prevents replacement.
outputMode = "overwrite";

% State dates: final calendar day of each month (inclusive).
startYear = 2003;
startMonth = 12;
endYear = 2024;
endMonth = 12;

% Parallel processing settings.
useParallel = true;
numWorkers = 15;

% Check the raw latitude/longitude order of each selected daily file.
checkAllStateCoordinates = true;

% Diagnostic map.
makeDiagnosticPlot = true;
plotDate = 20031231;

%% =========================================================================
% 3. MAIN EXECUTION
% =========================================================================

%% 3.1 Check paths and parameters
if ~isfolder(folderPath)
    error('Input directory not found: %s', folderPath);
end
if ~ismember(outputMode, ["overwrite", "stop"])
    error('outputMode must be "overwrite" or "stop".');
end
if useParallel && (~isscalar(numWorkers) || ~isfinite(numWorkers) || ...
        numWorkers < 1 || numWorkers ~= floor(numWorkers))
    error('numWorkers must be a positive integer.');
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

%% 3.2 Define monthly-end state dates
ymState = makeYMList(startYear, startMonth, endYear, endMonth);
nState = numel(ymState);
yearVecState = floor(ymState / 100);
monthVecState = mod(ymState, 100);
dayVecState = zeros(nState, 1);
for t = 1:nState
    dayVecState(t) = getMonthDays(yearVecState(t), monthVecState(t));
end
yyyymmddState = yearVecState * 10000 + monthVecState * 100 + dayVecState;
dateState = datetime(yearVecState, monthVecState, dayVecState, ...
    'Format', 'yyyy-MM-dd');

fprintf('GLDAS-DA V2.2 monthly-end GWS extraction\n');
fprintf('State dates: %d (%s to %s)\n', ...
    nState, char(dateState(1)), char(dateState(end)));

%% 3.3 Locate daily NetCDF files
ncFiles = cell(nState, 1);
missingDates = strings(0, 1);
for t = 1:nState
    thisDate = yyyymmddState(t);
    yearFolder = fullfile(folderPath, sprintf('%04d', yearVecState(t)));
    fileBase = sprintf('GLDAS_CLSM025_DA1_D.A%08d.022', thisDate);
    candidate1 = fullfile(yearFolder, [fileBase, '.nc4']);
    candidate2 = fullfile(yearFolder, fileBase);

    if isfile(candidate1)
        ncFiles{t} = candidate1;
    elseif isfile(candidate2)
        ncFiles{t} = candidate2;
    else
        % Select only NetCDF data files, not XML or checksum sidecars.
        listing = dir(fullfile(yearFolder, [fileBase, '*']));
        listing = listing(~[listing.isdir]);
        keep = ~cellfun('isempty', regexp({listing.name}, ...
            '^GLDAS_CLSM025_DA1_D\.A\d{8}\.022(\.nc4)?$', 'once'));
        listing = listing(keep);
        if numel(listing) == 1
            ncFiles{t} = fullfile(yearFolder, listing(1).name);
        elseif isempty(listing)
            missingDates(end + 1, 1) = string(sprintf('%08d', thisDate));
        else
            error('Multiple daily data files for %08d in %s', thisDate, yearFolder);
        end
    end
end
if ~isempty(missingDates)
    fprintf('Missing dates:\n%s\n', strjoin(missingDates, ', '));
    error('Required GLDAS-DA daily files are missing.');
end
fprintf('All %d daily files were found.\n', nState);

%% 3.4 Read coordinates and check the source grid
firstFile = ncFiles{1};
lat0 = double(ncread(firstFile, 'lat'));
lon0 = double(ncread(firstFile, 'lon'));
lon0(lon0 > 180) = lon0(lon0 > 180) - 360;

[lat, latOrder] = sort(lat0(:), 'ascend');
[lon, lonOrder] = sort(lon0(:), 'ascend');
nLat = numel(lat);
nLon = numel(lon);
validateNativeGrid(lat, lon);
[Lon2D, Lat2D] = meshgrid(lon, lat);
checkVariablesExist(firstFile, {'GWS_tavg'});

if checkAllStateCoordinates
    fprintf('Checking coordinates for %d monthly-end files...\n', nState);
    for t = 2:nState
        checkFileCoordinates(ncFiles{t}, lat0(:), lon0(:));
    end
    fprintf('All monthly-end coordinate grids are consistent.\n');
end
fprintf('Native grid: %d latitude x %d longitude.\n', nLat, nLon);

%% 3.5 Extract daily GWS_tavg for each selected date
GWS_mm = NaN(nLat, nLon, nState);
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
            'Parallel pool unavailable; continuing serially. Reason: %s', ME.message);
        useParallel = false;
    end
end

fprintf('Extracting GWS_tavg...\n');
runTimer = tic;
if useParallel
    parfor t = 1:nState
        GWS_mm(:, :, t) = readClsm2D( ...
            ncFiles{t}, 'GWS_tavg', nLat, nLon, latOrder, lonOrder);
    end
else
    for t = 1:nState
        fprintf('Processing %3d/%3d: %08d\n', t, nState, yyyymmddState(t));
        GWS_mm(:, :, t) = readClsm2D( ...
            ncFiles{t}, 'GWS_tavg', nLat, nLon, latOrder, lonOrder);
    end
end
fprintf('Extraction completed in %.2f minutes.\n', toc(runTimer) / 60);

%% 3.6 Check output
expectedSize = [nLat, nLon, nState];
checkOutputSize(GWS_mm, expectedSize, 'GWS_mm');
printVarSummary(GWS_mm, 'GWS_mm');

%% 3.7 Record metadata
Meta = struct();
Meta.sourceProduct = 'GLDAS_CLSM025_DA1_D_2.2';
Meta.sourceDOI = '10.5067/TXBMLX370XX8';
Meta.period = sprintf('%s to %s', char(dateState(1)), char(dateState(end)));
Meta.variable = 'GWS_tavg';
Meta.unit = 'mm';
Meta.arrayConvention = '[lat, lon, time], lat south-to-north, lon west-to-east';
Meta.temporalDefinition = ...
    'Daily mean GWS_tavg from the final calendar day of each month';
Meta.monthlyChangeDefinition = 'dGWS(t) = GWS_monthend(t) - GWS_monthend(t-1)';
Meta.referenceMonth = '2003-12';
Meta.originalResolutionDegrees = 0.25;
Meta.spatialUpscalingApplied = false;
Meta.coordinateConsistencyChecked = checkAllStateCoordinates;
Meta.nSourceDates = nState;
Meta.inputFilePattern = 'GLDAS_CLSM025_DA1_D.AYYYYMMDD.022.nc4';
Meta.documentation = 'NASA README_GLDAS2.pdf, Sections 3.1 and 3.5';

%% 3.8 Save output
fprintf('Saving MAT-file:\n  %s\n', outFile);
save(outFile, ...
    'GWS_mm', ...
    'lat', 'lon', 'Lat2D', 'Lon2D', ...
    'dateState', 'ymState', 'yyyymmddState', ...
    'yearVecState', 'monthVecState', 'dayVecState', ...
    'Meta', ...
    '-v7.3');
fprintf('Saved successfully.\n');

%% 3.9 Diagnostic map
if makeDiagnosticPlot
    plotIndex = find(yyyymmddState == plotDate, 1);
    if isempty(plotIndex)
        warning('SKLGB:PlotDateMissing', ...
            'Diagnostic map skipped: %08d is not in the extracted range.', plotDate);
    else
        gwsMap = GWS_mm(:, :, plotIndex);
        figure('Name', 'GLDAS-DA 0.25-degree GWS check');
        p = imagesc(lon, lat, double(gwsMap));
        set(p, 'AlphaData', double(isfinite(gwsMap)));
        set(gca, 'YDir', 'normal');
        axis equal tight;
        box on;
        grid on;
        colormap(jet);
        cb = colorbar;
        ylabel(cb, 'GWS (mm)');
        title(['GLDAS-DA CLSM GWS, ', char(dateState(plotIndex)), ...
            ', native 0.25^{\circ}']);
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
        error('End month cannot precede start month.');
    end
    d = d1:calmonths(1):d2;
    ym = double(year(d(:)) * 100 + month(d(:)));
end

function nDays = getMonthDays(y, m)
    y = double(y);
    m = double(m);
    if ~isscalar(y) || ~isscalar(m) || ~isfinite(y) || ~isfinite(m) || ...
            y ~= floor(y) || m ~= floor(m) || m < 1 || m > 12
        error('Year and month must be valid integer scalars.');
    end
    d1 = datetime(y, m, 1);
    d2 = d1 + calmonths(1);
    nDays = double(days(d2 - d1));
    if ~isscalar(nDays) || ~isfinite(nDays) || nDays < 28 || nDays > 31
        error('Unexpected days-in-month value for %04d-%02d.', y, m);
    end
end

function validateNativeGrid(lat, lon)
    if isempty(lat) || isempty(lon) || ...
            any(~isfinite(lat)) || any(~isfinite(lon)) || ...
            any(diff(lat) <= 0) || any(diff(lon) <= 0)
        error('Native latitude/longitude coordinates are invalid.');
    end
    tol = 1e-5;
    if numel(lat) ~= 600 || numel(lon) ~= 1440
        error('Expected the full 600 x 1440 GLDAS-DA source grid.');
    end
    if any(abs(diff(lat) - 0.25) > tol) || ...
            any(abs(diff(lon) - 0.25) > tol)
        error('Expected a regular 0.25-degree source grid.');
    end
    if min(lat) < -90 || max(lat) > 90 || min(lon) < -180 || max(lon) > 180
        error('Latitude or longitude is outside the expected range.');
    end
end

function checkFileCoordinates(ncFile, expectedLatRaw, expectedLonRaw)
    thisLat = double(ncread(ncFile, 'lat'));
    thisLon = double(ncread(ncFile, 'lon'));
    thisLon(thisLon > 180) = thisLon(thisLon > 180) - 360;
    thisLat = thisLat(:);
    thisLon = thisLon(:);
    if numel(thisLat) ~= numel(expectedLatRaw) || ...
            numel(thisLon) ~= numel(expectedLonRaw) || ...
            any(~isfinite(thisLat)) || any(~isfinite(thisLon)) || ...
            any(abs(thisLat - expectedLatRaw) > 1e-5) || ...
            any(abs(thisLon - expectedLonRaw) > 1e-5)
        error('Monthly-end NetCDF grid does not match the first file: %s', ncFile);
    end
end

function checkVariablesExist(ncFile, varNames)
    info = ncinfo(ncFile);
    existing = {info.Variables.Name};
    for i = 1:numel(varNames)
        if ~any(strcmp(existing, varNames{i}))
            error('Variable "%s" not found in file: %s', varNames{i}, ncFile);
        end
    end
end

function A = readClsm2D(ncFile, varName, nLat, nLon, latOrder, lonOrder)
    raw = squeeze(double(ncread(ncFile, varName)));
    fillVals = [];
    try
        att = double(ncreadatt(ncFile, varName, '_FillValue'));
        fillVals = [fillVals; att(:)];
    catch
    end
    try
        att = double(ncreadatt(ncFile, varName, 'missing_value'));
        fillVals = [fillVals; att(:)];
    catch
    end
    for k = 1:numel(fillVals)
        raw(raw == fillVals(k)) = NaN;
    end
    raw(~isfinite(raw) | abs(raw) > 1e19) = NaN;

    if isequal(size(raw), [nLon, nLat])
        A = raw';
    elseif isequal(size(raw), [nLat, nLon])
        A = raw;
    else
        error('Unexpected size for "%s" in "%s": [%s].', ...
            varName, ncFile, num2str(size(raw)));
    end
    A = A(latOrder, lonOrder);
end

function checkOutputSize(X, expectedSize, varName)
    if ~isequal(size(X), expectedSize) || ~isa(X, 'double')
        error('Unexpected output size or class for %s.', varName);
    end
end

function printVarSummary(X, varName)
    fprintf('%-10s | class: %-6s | size: %-18s | NaN frac: %.4f | min: %.4f | mean: %.4f | max: %.4f\n', ...
        varName, class(X), mat2str(size(X)), mean(isnan(X(:))), ...
        min(X(:), [], 'omitnan'), mean(X(:), 'omitnan'), ...
        max(X(:), [], 'omitnan'));
end
