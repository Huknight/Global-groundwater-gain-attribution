%% STEP 1 - Extract GLDAS-Noah 0.25-degree monthly variables
% =========================================================================
% 1. CODE DESCRIPTION
% =========================================================================
% Purpose
%   Extract five monthly hydroclimatic fields from NASA
%   GLDAS-Noah Version 2.1 0.25-degree monthly NetCDF files (Main Production).
%
% Data period
%   December 2003 to December 2024, inclusive (253 monthly records).
%   December 2003 is retained for the subsequent soil-moisture difference.
%
% Input product
%   GLDAS_NOAH025_M_2.1, filenames GLDAS_NOAH025_M.AYYYYMM.021.nc4
%   NASA product DOI: 10.5067/SXAVCZFAQLNO
%
% Output variables and conversions (GLDAS Version 2 README, Sections 3.3/3.5)
%   Prec_mm   = Rainf_f_tavg * 86400 * days_in_month       [mm/month]
%   Evap_mm   = Evap_tavg    * 86400 * days_in_month       [mm/month]
%   Runoff_mm = (Qs_acc + Qsb_acc) * 8 * days_in_month     [mm/month]
%   SM2m_mm   = sum of four Noah soil-moisture layers      [mm]
%   SurfT_C   = AvgSurfT_inst - 273.15                     [deg C]
%
% Notes on physical meaning
%   * Rainf_f_tavg is TOTAL precipitation (liquid + frozen).
%   * Qs_acc and Qsb_acc are mean 3-hour ACCUMULATIONS in the monthly file.
%   * Soil moisture is the total 0-200 cm water-equivalent storage, not a
%     root-zone fraction or volumetric moisture content.
%   * The monthly AvgSurfT_inst field is an average of instantaneous
%     surface-skin temperature samples, not an end-of-month temperature.
%   * 1 kg/m^2 of water is numerically equivalent to 1 mm water depth.
%
% Output MAT contents
%   Prec_mm, Evap_mm, Runoff_mm, SM2m_mm, SurfT_C : [lat, lon, time]
%   lat, lon                                  : ascending coordinate vectors
%   Lat2D, Lon2D                              : coordinate grids
%   ym, yearVec, monthVec                     : month identifiers
%   Meta                                      : processing documentation
%   MAT-file format                           : -v7.3
%
% Reference
%   NASA GES DISC, README_GLDAS2.pdf, particularly Sections 3.3 and 3.5.
%

%% =========================================================================
% 2. CENTRALIZED PARAMETERS
% =========================================================================
clearvars;
clc;

% Input NetCDF directory.
folderPath = ...
    'D:\A1_科研数据库\01_原始数据\01_GLDAS\GLDAS2.1 Noah Monthly 025D 2000-01-to-2025-12';

% Output directory.
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised\s1_master_database';
outFolder = fullfile(projectRoot, '01_GLDAS_extraction');
outFileName = 'GLDAS_NOAH025_native5_200312_202412.mat';

% Output file handling: "overwrite" replaces the file; "stop" prevents replacement.
outputMode = "overwrite";

% Data period (inclusive).
startYear = 2003;
startMonth = 12;
endYear = 2024;
endMonth = 12;

% Parallel processing settings.
useParallel = true;
numWorkers = 15;

% Check latitude and longitude consistency across monthly files.
checkAllMonthlyCoordinates = true;

% Diagnostic precipitation map.
makeDiagnosticPlot = true;
plotYM = 200312;

%% =========================================================================
% 3. MAIN EXECUTION
% =========================================================================

%% 3.1 Validate settings and prepare paths
if ~isfolder(folderPath)
    error('Input directory not found: %s', folderPath);
end

if ~ismember(outputMode, ["overwrite", "stop"])
    error('outputMode must be "overwrite" or "stop".');
end

if useParallel && (~isscalar(numWorkers) || numWorkers < 1 || ...
        numWorkers ~= floor(numWorkers))
    error('numWorkers must be a positive integer.');
end

if ~exist(outFolder, 'dir')
    [madeFolder, message] = mkdir(outFolder);
    if ~madeFolder
        error('Cannot create output directory: %s\n%s', outFolder, message);
    end
end
outFile = fullfile(outFolder, outFileName);

if exist(outFile, 'file') && outputMode == "stop"
    error('Output file already exists: %s', outFile);
end

%% 3.2 Define the exact monthly sequence
ym = makeYMList(startYear, startMonth, endYear, endMonth);
nTime = numel(ym);
yearVec = floor(ym / 100);
monthVec = mod(ym, 100);

fprintf('GLDAS-Noah V2.1 monthly extraction\n');
fprintf('Months: %d (%04d-%02d to %04d-%02d)\n', ...
    nTime, yearVec(1), monthVec(1), yearVec(end), monthVec(end));

%% 3.3 Identify and validate source files
files = dir(fullfile(folderPath, 'GLDAS_NOAH025_M.A*.021.nc4'));
if isempty(files)
    error('No GLDAS_NOAH025_M.A*.021.nc4 files found in: %s', folderPath);
end

fileNames = {files.name}';
fileYM = nan(numel(fileNames), 1);
for k = 1:numel(fileNames)
    tok = regexp(fileNames{k}, ...
        '^GLDAS_NOAH025_M\.A(\d{6})\.021\.nc4$', 'tokens', 'once');
    if ~isempty(tok)
        fileYM(k) = str2double(tok{1});
    end
end

% Keep files matching the GLDAS-Noah monthly naming convention.
validNames = ~isnan(fileYM);
files = files(validNames);
fileYM = fileYM(validNames);

[fileYM, sortOrder] = sort(fileYM);
files = files(sortOrder);

if numel(unique(fileYM)) ~= numel(fileYM)
    error('Duplicate YYYYMM records among the source NetCDF files.');
end

missingYM = setdiff(ym(:), fileYM(:));
if ~isempty(missingYM)
    fprintf('Missing YYYYMM: %s\n', num2str(missingYM(:)'));
    error('Required GLDAS-Noah monthly files are missing.');
end

[foundEveryMonth, selectedIdx] = ismember(ym(:), fileYM(:));
if ~all(foundEveryMonth)
    error('Internal error while mapping the required months to files.');
end

ncFiles = cell(nTime, 1);
for t = 1:nTime
    ncFiles{t} = fullfile(folderPath, files(selectedIdx(t)).name);
end
fprintf('All %d required monthly NetCDF files are present.\n', nTime);

%% 3.4 Read coordinate vectors and check variable availability
firstFile = ncFiles{1};
lat0 = double(ncread(firstFile, 'lat'));
lon0 = double(ncread(firstFile, 'lon'));

% Convert longitude to [-180, 180] if needed.
lon0(lon0 > 180) = lon0(lon0 > 180) - 360;

% Output orientation: latitude south-to-north, longitude west-to-east.
[lat, latOrder] = sort(lat0(:), 'ascend');
[lon, lonOrder] = sort(lon0(:), 'ascend');
nLat = numel(lat);
nLon = numel(lon);

validateNativeGrid(lat, lon);
[Lon2D, Lat2D] = meshgrid(lon, lat);
fprintf('Native spatial dimensions: %d latitude x %d longitude.\n', nLat, nLon);

requiredVars = { ...
    'Rainf_f_tavg', ...
    'Evap_tavg', ...
    'Qs_acc', ...
    'Qsb_acc', ...
    'SoilMoi0_10cm_inst', ...
    'SoilMoi10_40cm_inst', ...
    'SoilMoi40_100cm_inst', ...
    'SoilMoi100_200cm_inst', ...
    'AvgSurfT_inst'};
checkVariablesExist(firstFile, requiredVars);

% Verify monthly spatial coordinates before reading data fields.
if checkAllMonthlyCoordinates
    fprintf('Checking monthly latitude/longitude consistency...\n');
    for t = 2:nTime
        checkFileCoordinates(ncFiles{t}, lat0(:), lon0(:));
    end
    fprintf('All %d monthly coordinate grids are consistent.\n', nTime);
end

%% 3.5 Allocate output arrays
Prec_mm   = NaN(nLat, nLon, nTime);
Evap_mm   = NaN(nLat, nLon, nTime);
Runoff_mm = NaN(nLat, nLon, nTime);
SM2m_mm   = NaN(nLat, nLon, nTime);
SurfT_C   = NaN(nLat, nLon, nTime);

%% 3.6 Start parallel processing if requested and available
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
            ['Parallel pool unavailable; continuing serially. ' ...
             'Reason: %s'], ME.message);
        useParallel = false;
    end
end

%% 3.7 Extract and convert the monthly fields
fprintf('Beginning extraction (%s mode).\n', ternaryMode(useParallel));
tic;

if useParallel
    parfor t = 1:nTime
        daysThisMonth = getMonthDays(yearVec(t), monthVec(t));
        [Prec_t, Evap_t, Runoff_t, SM2m_t, SurfT_t] = extractOneMonthNoah( ...
            ncFiles{t}, nLat, nLon, latOrder, lonOrder, daysThisMonth);

        Prec_mm(:, :, t)   = Prec_t;
        Evap_mm(:, :, t)   = Evap_t;
        Runoff_mm(:, :, t) = Runoff_t;
        SM2m_mm(:, :, t)   = SM2m_t;
        SurfT_C(:, :, t)   = SurfT_t;
    end
else
    for t = 1:nTime
        daysThisMonth = getMonthDays(yearVec(t), monthVec(t));
        fprintf('Processing %3d/%3d: %06d\n', t, nTime, ym(t));
        [Prec_t, Evap_t, Runoff_t, SM2m_t, SurfT_t] = extractOneMonthNoah( ...
            ncFiles{t}, nLat, nLon, latOrder, lonOrder, daysThisMonth);

        Prec_mm(:, :, t)   = Prec_t;
        Evap_mm(:, :, t)   = Evap_t;
        Runoff_mm(:, :, t) = Runoff_t;
        SM2m_mm(:, :, t)   = SM2m_t;
        SurfT_C(:, :, t)   = SurfT_t;
    end
end

elapsedTime = toc;
fprintf('Extraction completed in %.2f minutes.\n', elapsedTime / 60);

%% 3.8 Record metadata
Meta = struct();
Meta.sourceProduct = 'GLDAS_NOAH025_M_2.1';
Meta.period = sprintf('%04d-%02d to %04d-%02d', ...
    startYear, startMonth, endYear, endMonth);
Meta.arrayConvention = '[lat, lon, time], lat south-to-north, lon west-to-east';
Meta.Prec_mm = 'Rainf_f_tavg * 86400 * days_in_month; unit: mm/month';
Meta.Evap_mm = 'Evap_tavg * 86400 * days_in_month; unit: mm/month';
Meta.Runoff_mm = '(Qs_acc + Qsb_acc) * 8 * days_in_month; unit: mm/month';
Meta.SM2m_mm = ['SoilMoi0_10cm_inst + SoilMoi10_40cm_inst + ' ...
    'SoilMoi40_100cm_inst + SoilMoi100_200cm_inst; unit: mm'];
Meta.SurfT_C = 'AvgSurfT_inst - 273.15; unit: degree C';
Meta.note1 = 'Qsm_acc is not included in Runoff_mm.';

% Source and processing metadata.
Meta.sourceDOI = '10.5067/SXAVCZFAQLNO';
Meta.inputFilePattern = 'GLDAS_NOAH025_M.AYYYYMM.021.nc4';
Meta.originalResolutionDegrees = 0.25;
Meta.outputResolutionDegrees = 0.25;
Meta.nSourceMonths = nTime;
Meta.coordinateConsistencyChecked = checkAllMonthlyCoordinates;
Meta.spatialUpscalingApplied = false;
Meta.scriptName = 'step1_extract_noah_data.m';
Meta.documentation = 'NASA README_GLDAS2.pdf, Sections 3.3 and 3.5';

%% 3.9 Save MAT output
fprintf('Saving MAT-file:\n  %s\n', outFile);
% Save the completed extraction output.
save(outFile, ...
    'Prec_mm', 'Evap_mm', 'Runoff_mm', 'SM2m_mm', 'SurfT_C', ...
    'lat', 'lon', 'Lat2D', 'Lon2D', ...
    'ym', 'yearVec', 'monthVec', 'Meta', ...
    '-v7.3');
fprintf('Saved successfully.\n');

%% 3.10 Plot precipitation for spatial coverage inspection
if makeDiagnosticPlot
    plotIndex = find(ym == plotYM, 1);
    if isempty(plotIndex)
        warning('SKLGB:PlotMonthMissing', ...
            'Diagnostic plot skipped: %06d not in the extracted range.', plotYM);
    else
        rainMap = Prec_mm(:, :, plotIndex);
        figure('Name', 'GLDAS-Noah 0.25-degree precipitation check');
        % Display the complete native grid.
        p = imagesc(lon, lat, double(rainMap));
        set(p, 'AlphaData', double(isfinite(rainMap)));
        set(gca, 'YDir', 'normal');
        axis equal tight;
        box on;
        grid on;
        colormap(jet);
        cb = colorbar;
        ylabel(cb, 'Precipitation (mm/month)');
        title(sprintf('GLDAS-Noah total precipitation, %04d-%02d, native 0.25^{\\circ}', ...
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
        error('End month cannot precede start month.');
    end
    d = d1:calmonths(1):d2;
    ym = double(year(d(:)) * 100 + month(d(:)));
end

function nDays = getMonthDays(y, m)
    y = double(y);
    m = double(m);
    if ~isscalar(y) || ~isscalar(m) || ~isfinite(y) || ...
            ~isfinite(m) || y ~= floor(y) || m ~= floor(m) || ...
            m < 1 || m > 12
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
    % Detect unexpected resampling, invalid coordinates, or misordered axes.
    if isempty(lat) || isempty(lon) || ...
            any(~isfinite(lat)) || any(~isfinite(lon)) || ...
            any(diff(lat) <= 0) || any(diff(lon) <= 0)
        error('Native latitude/longitude coordinates are invalid.');
    end
    tol = 1e-5;
    % The full GLDAS-Noah 0.25-degree grid has 600 x 1440 cells.
    if numel(lat) ~= 600 || numel(lon) ~= 1440
        error('Expected the full 600 x 1440 GLDAS-Noah source grid.');
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
    % Compare coordinates in source order to detect axis reversal.
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
        error('Monthly NetCDF grid does not match the first month: %s', ncFile);
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

function modeName = ternaryMode(isParallel)
    if isParallel
        modeName = 'parallel';
    else
        modeName = 'serial';
    end
end

function [Prec_t, Evap_t, Runoff_t, SM2m_t, SurfT_t] = extractOneMonthNoah( ...
    ncFile, nLat, nLon, latOrder, lonOrder, daysThisMonth)

    rainRate = readNoah2D(ncFile, 'Rainf_f_tavg', nLat, nLon, latOrder, lonOrder);
    evapRate = readNoah2D(ncFile, 'Evap_tavg', nLat, nLon, latOrder, lonOrder);

    qsAcc = readNoah2D(ncFile, 'Qs_acc', nLat, nLon, latOrder, lonOrder);
    qsbAcc = readNoah2D(ncFile, 'Qsb_acc', nLat, nLon, latOrder, lonOrder);

    sm1 = readNoah2D(ncFile, 'SoilMoi0_10cm_inst', nLat, nLon, latOrder, lonOrder);
    sm2 = readNoah2D(ncFile, 'SoilMoi10_40cm_inst', nLat, nLon, latOrder, lonOrder);
    sm3 = readNoah2D(ncFile, 'SoilMoi40_100cm_inst', nLat, nLon, latOrder, lonOrder);
    sm4 = readNoah2D(ncFile, 'SoilMoi100_200cm_inst', nLat, nLon, latOrder, lonOrder);

    surfT_K = readNoah2D(ncFile, 'AvgSurfT_inst', nLat, nLon, latOrder, lonOrder);

    % Monthly water-depth and temperature conversions.
    Prec_t = rainRate * 86400 * daysThisMonth;
    Evap_t = evapRate * 86400 * daysThisMonth;
    Runoff_t = (qsAcc + qsbAcc) * 8 * daysThisMonth;
    SM2m_t = sm1 + sm2 + sm3 + sm4;
    SurfT_t = surfT_K - 273.15;
end

function A = readNoah2D(ncFile, varName, nLat, nLon, latOrder, lonOrder)
    raw = double(ncread(ncFile, varName));
    raw = squeeze(raw);

    % ncread handles _FillValue and scale_factor/add_offset. Check for
    % remaining missing-value sentinels and invalid values.
    fillVals = [];
    try
        att = double(ncreadatt(ncFile, varName, '_FillValue'));
        fillVals = [fillVals; att(:)];
    catch
        % Attribute is optional.
    end
    try
        att = double(ncreadatt(ncFile, varName, 'missing_value'));
        fillVals = [fillVals; att(:)];
    catch
        % Attribute is optional.
    end
    for k = 1:numel(fillVals)
        raw(raw == fillVals(k)) = NaN;
    end
    raw(abs(raw) > 1e19) = NaN;

    % Arrange all fields as [lat, lon].
    if isequal(size(raw), [nLon, nLat])
        A = raw';
    elseif isequal(size(raw), [nLat, nLon])
        A = raw;
    else
        error(['Unexpected size of "%s" in "%s": [%s]. ' ...
            'Expected [%d %d] or [%d %d].'], ...
            varName, ncFile, num2str(size(raw)), nLon, nLat, nLat, nLon);
    end
    A = A(latOrder, lonOrder);
end
