%% koppen_climate_attribution.m
% Regional gain-attribution analysis by Koppen-Geiger climate zones.
%
% INPUT 1: s3_uc_analysis1_results.mat (original run, 2022-2024)
% INPUT 2: Koppen_1991_2020/koppen_geiger_1p0.tif (Beck et al. 2023)
%
% This self-contained script:
%   (1) matches 17,284 valid GLDAS grid centres to 1-degree climate classes;
%   (2) RECONSTRUCTS the primary three-year-mean-NRMSE attribution using
%       the rules of s2_gain_attribution2.m;
%   (3) checks the reconstructed primary class counts against the manuscript;
%   (4) computes the same statistics for each year (2022, 2023, 2024),
%       checking annual maps against s3_uc_analysis1_results.mat;
%   (5) summarizes six gain categories, nine winners and component shares
%       within each of five major Koppen-Geiger zones (A-E);
%   (6) saves MAT, CSV, optional XLSX, and PNG figures.
%
% The five climate zones are A=1-3, B=4-7, C=8-16, D=17-28, E=29-30.
% Within-hydroclimatic-gain component shares use SUMMED NRMSE reductions,
% not averages of individual grid-cell percentages.
%
% No retraining, interpolation, or use of target-period predictor data.
% Requires MATLAB R2018b or later for string/table/export functions.
% readgeoraster is used when available; imread is a fallback for this
% verified 180-by-360 GeoTIFF (standard global geographic grid).
%
% Data source: Beck et al. (2023), Scientific Data 10:724,
% DOI: 10.1038/s41597-023-02549-6; corrected Figshare release (2026).

clearvars; clc;

%% ========================= USER SETTINGS ===============================
% If files are not next to the script, leave blank for file picker.
% Alternatively provide absolute file paths directly in these variables.
matFile = '';
tifFile = '';

% Portable default output directory alongside this script.
% All CSV, XLSX, MAT and PNG outputs are saved here.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
outDir = fullfile(scriptDir, 'outputs');

saveFigures = true;
writeExcel = true;       % CSV files are always saved even if XLSX fails.

%% ========================= FIND INPUT FILES ============================

if isempty(matFile)
    candidates = { ...
        fullfile(scriptDir, 's3_uc_analysis1_results.mat'), ...
        fullfile(scriptDir, 's3_uc_analysis1_results', 's3_uc_analysis1_results.mat'), ...
        fullfile(scriptDir, 's3_darft', 's1_analysis', ...
                 's3_uc_analysis1_results', 's3_uc_analysis1_results.mat')};
    matFile = firstExisting(candidates);
end
if isempty(matFile)
    [name, folder] = uigetfile('*.mat', ...
        'Select s3_uc_analysis1_results.mat');
    if isequal(name, 0), error('No MAT input selected.'); end
    matFile = fullfile(folder, name);
end
if ~isfile(matFile), error('MAT file not found: %s', matFile); end

if isempty(tifFile)
    candidates = { ...
        fullfile(scriptDir, 'koppen_geiger_1p0.tif'), ...
        fullfile(scriptDir, 'Koppen_1991_2020', 'koppen_geiger_1p0.tif'), ...
        fullfile(fileparts(matFile), 'koppen_geiger_1p0.tif'), ...
        fullfile(fileparts(matFile), 'Koppen_1991_2020', 'koppen_geiger_1p0.tif')};
    tifFile = firstExisting(candidates);
end
if isempty(tifFile)
    [name, folder] = uigetfile({'*.tif;*.tiff', 'GeoTIFF files'}, ...
        'Select 1991_2020/koppen_geiger_1p0.tif');
    if isequal(name, 0), error('No GeoTIFF input selected.'); end
    tifFile = fullfile(folder, name);
end
if ~isfile(tifFile), error('GeoTIFF file not found: %s', tifFile); end

if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('MAT: %s\nTIFF: %s\nOUTPUT: %s\n\n',matFile,tifFile,outDir);

%% ========================= LOAD RESULTS ================================
S = load(matFile, 'lat', 'lon', 'commonTemporalMask', ...
    'annualNRMSEStack', 'annualResult', 'testYears', 'modelNames');
required = {'lat', 'lon', 'commonTemporalMask', ...
    'annualNRMSEStack', 'annualResult', 'testYears', 'modelNames'};
for i = 1:numel(required)
    if ~isfield(S, required{i})
        error('Missing required variable in MAT: %s', required{i});
    end
end

lat = double(S.lat(:));
lon = double(S.lon(:));
validMask = logical(S.commonTemporalMask);
annualNRMSEStack = double(S.annualNRMSEStack);
years = double(S.testYears(:))';
modelNames = string(S.modelNames(:));

nLat = numel(lat); nLon = numel(lon);
if ~isequal(size(validMask), [nLat, nLon])
    error('Valid mask dimensions disagree with lat/lon vectors.');
end
if ~isequal(size(annualNRMSEStack), [nLat, nLon, 9, 3])
    error('Expected annualNRMSEStack to be [nLat,nLon,9,3]; actual: %s', ...
        mat2str(size(annualNRMSEStack)));
end
if numel(years) ~= 3 || ~isequal(years, [2022, 2023, 2024])
    error('Expected testYears [2022 2023 2024]; actual: %s',mat2str(years));
end
if numel(modelNames) ~= 9
    error('Expected 9 model names.');
end
if nnz(validMask) ~= 17284
    error('Expected 17284 valid cells, obtained %d.', nnz(validMask));
end
if any(~isfinite(annualNRMSEStack(repmat(validMask,1,1,9,3))))
    error('At least one valid grid cell has non-finite annual NRMSE.');
end

fprintf('Checked: %d valid cells, %d models, test years %s.\n', ...
    nnz(validMask), numel(modelNames), mat2str(years));

%% ========================= GEOREFERENCED MATCH =========================
if exist('readgeoraster', 'file') == 2
    [K, R] = readgeoraster(tifFile);
    if isprop(R, 'LatitudeLimits') && isprop(R, 'LongitudeLimits')
        if max(abs(double(R.LatitudeLimits) - [-90, 90])) > 1e-6 || ...
           max(abs(double(R.LongitudeLimits) - [-180, 180])) > 1e-6
            error('GeoTIFF spatial bounds are not global [-180,180]/[-90,90].');
        end
    end
else
    % This fallback is valid ONLY for the named, verified 1-degree map.
    [~, tifName, tifExt] = fileparts(tifFile);
    if ~strcmpi([tifName, tifExt], 'koppen_geiger_1p0.tif')
        warning('Verify the selected TIFF is the 1-degree 1991-2020 map.');
    end
    K = imread(tifFile);
end
if ndims(K) ~= 2 || ~isequal(size(K), [180, 360])
    error('Expected 180x360 single-band Koppen raster, got %s.',mat2str(size(K)));
end
K = double(K);

% The 1991-2020 1-degree TIFF is north-to-south, west-to-east.
% Original GLDAS latitude is ascending (-59.5 to 89.5), so DO NOT align
% arrays directly. Match each GLDAS centre using its actual coordinates.
lon180 = mod(lon + 180, 360) - 180;
if any(abs(lat - (floor(lat) + 0.5)) > 1e-7) || ...
   any(abs(lon180 - (floor(lon180) + 0.5)) > 1e-7)
    error('GLDAS centres are not at half-degree grid centres.');
end
if any(lat <= -90 | lat >= 90) || any(lon180 < -180 | lon180 >= 180)
    error('Coordinates outside the reference raster.');
end

rowIdx = floor(90 - lat) + 1;
colIdx = floor(lon180 + 180) + 1;
if any(rowIdx < 1 | rowIdx > 180) || any(colIdx < 1 | colIdx > 360)
    error('Grid-to-raster indices out of bounds.');
end
kgFineMap = K(rowIdx, colIdx);    % [150 x 360], same as validMask
kgFineMap(~validMask) = NaN;

kgMajorMap = nan(nLat, nLon);
kgMajorMap(kgFineMap >= 1  & kgFineMap <= 3 ) = 1; % A
kgMajorMap(kgFineMap >= 4  & kgFineMap <= 7 ) = 2; % B
kgMajorMap(kgFineMap >= 8  & kgFineMap <= 16) = 3; % C
kgMajorMap(kgFineMap >= 17 & kgFineMap <= 28) = 4; % D
kgMajorMap(kgFineMap >= 29 & kgFineMap <= 30) = 5; % E

nClimateValid = nnz(validMask & isfinite(kgMajorMap));
if nClimateValid ~= nnz(validMask)
    error('Only %d/%d valid GLDAS cells obtained climate labels. Check input.', ...
        nClimateValid, nnz(validMask));
end
climateLetters = ["A", "B", "C", "D", "E"];
climateNames = ["Tropical", "Arid", "Temperate", "Cold", "Polar"];
climateCounts = zeros(5,1);
for z=1:5
    climateCounts(z) = nnz(validMask & kgMajorMap == z);
end
fprintf('Climate matching COMPLETE (A-E): %s\n\n',mat2str(climateCounts'));

%% ========================= PRIMARY + ANNUAL ============================
periodLabels = ["ThreeYearMean", string(years)];
results = struct([]);
validationRecords = table();

for p=1:4
    if p==1
        % Manuscript primary attribution: mean model NRMSE across THREE
        % 12-month test years, THEN select winner and decompose gains.
        modelNRMSE = mean(annualNRMSEStack, 4);
    else
        modelNRMSE = annualNRMSEStack(:,:,:,p-1);
    end

    R = deriveAttribution(modelNRMSE, validMask);
    % MATLAB does not allow indexed assignment of a populated struct into
    % struct([]) (which has no fields). Initialise from the first result;
    % subsequent results have the same fields and can be indexed normally.
    if p == 1
        results = R;
    else
        results(p) = R; %#ok<SAGROW>
    end

    % Reproduce the manuscript main classification exactly.
    if p==1
        refCounts = [1406, 4343, 4503, 2083, 3008, 1941];
        thisCounts = arrayfun(@(c)nnz(R.class6Map(validMask) == c),1:6);
        if ~isequal(thisCounts,refCounts)
            error(['Primary classification does not match Section 3.2. ', ...
                'Expected %s; obtained %s.'],mat2str(refCounts),mat2str(thisCounts));
        end
    else
        % Independent check with saved annual results from s3 script.
        saved = S.annualResult(p-1);
        cMismatch = nnz(validMask & R.gainClassMap ~= saved.gainClassMap);
        wMismatch = nnz(validMask & R.bestAll9Idx ~= saved.bestAll9Idx);
        if cMismatch>0 || wMismatch>0
            error(['Year %d reconstructed maps do not match stored annual ', ...
                'attribution: %d category, %d winner mismatches.'], ...
                years(p-1),cMismatch,wMismatch);
        end
    end

    vc = table(periodLabels(p), R.nValid, R.maxReconstructionError, ...
        'VariableNames', {'Period','NValid','MaxAbsGainReconstructionError'});
    validationRecords = [validationRecords; vc]; %#ok<AGROW>

    fprintf('%-15s: categories 1-6 = %s, reconstruction error %.2g\n',...
        periodLabels(p),mat2str(arrayfun(@(c)nnz(validMask & R.class6Map==c),1:6)),...
        R.maxReconstructionError);
end

%% ========================= SUMMARIZE BY CLIMATE =========================
zoneTable = table();
winnerTable = table();
categoryTable = table();

sixLabels = ["No gain", "Linear gain", "Nonlinear gain", ...
             "Hydroclimatic gain", "Mixed-L", "Mixed-N"];

for p=1:4
    R = results(p);
    for z=0:5
        if z==0
            zMask = validMask;
            zLetter = "GLOBAL";
            zName = "All valid grid cells";
        else
            zMask = validMask & kgMajorMap==z;
            zLetter = climateLetters(z);
            zName = climateNames(z);
        end
        n = nnz(zMask);
        if n==0, error('No cells in zone %s.',zLetter); end
        cCounts = arrayfun(@(c)nnz(zMask & R.class6Map==c),1:6);
        cPct = 100*cCounts/n;
        nHydro = cCounts(4);
        hydroMask = zMask & R.class6Map==4;
        totalHydroGain = sum(R.G_total(hydroMask));
        if nHydro>0 && totalHydroGain>0
            hHistorical = 100*sum(R.G_linear(hydroMask) + ...
                                   R.G_nonlinear(hydroMask))/totalHydroGain;
            hHydroclimatic = 100*sum(R.G_hydro(hydroMask))/totalHydroGain;
        else
            hHistorical = NaN;
            hHydroclimatic = NaN;
        end

        zr = table(periodLabels(p), zLetter, zName, n, ...
            100*n/nnz(validMask), ...
            cCounts(1), cCounts(2), cCounts(3), cCounts(4), ...
            cCounts(5), cCounts(6), ...
            cPct(1), cPct(2), cPct(3), cPct(4), cPct(5), cPct(6), ...
            cPct(5)+cPct(6), cPct(2)+cPct(3), ...
            nHydro, hHistorical, hHydroclimatic, ...
            mean(R.allFusion(zMask)), mean(R.AR_NRMSE(zMask)), ...
            mean(R.G_total(zMask)), ...
            'VariableNames', { ...
            'Period','ClimateZone','ClimateName','NCells','GlobalCellSharePct', ...
            'NoGainN','LinearN','NonlinearN','HydroclimaticN','MixedLN','MixedNN', ...
            'NoGainPct','LinearPct','NonlinearPct','HydroclimaticPct', ...
            'MixedLPct','MixedNPct','MixedTotalPct','StructuralCategoriesPct', ...
            'NHydroCells','HistoricalContributionPct','HydroclimaticContributionPct', ...
            'MeanWinnerNRMSE','MeanARNRMSE','MeanTotalGainNRMSE'});
        zoneTable = [zoneTable;zr]; %#ok<AGROW>

        for c=1:6
            cr = table(periodLabels(p),zLetter,sixLabels(c),cCounts(c),cPct(c), ...
                'VariableNames',{'Period','ClimateZone','GainCategory','NCells','PctWithinZone'});
            categoryTable = [categoryTable;cr]; %#ok<AGROW>
        end
        for m=1:9
            numW = nnz(zMask & R.bestAll9Idx==m);
            wr = table(periodLabels(p),zLetter,modelNames(m),numW,100*numW/n, ...
                'VariableNames',{'Period','ClimateZone','WinningModel','NCells','PctWithinZone'});
            winnerTable = [winnerTable;wr]; %#ok<AGROW>
        end
    end
end

% Verify zone counts and category/winner totals per period.
for p=1:4
    for z=1:5
        v = zoneTable.Period==periodLabels(p) & zoneTable.ClimateZone==climateLetters(z);
        if nnz(v)~=1 || zoneTable.NCells(v)~=climateCounts(z)
            error('Zone count consistency check failed: %s %s',periodLabels(p),climateLetters(z));
        end
    end
    v = zoneTable.Period==periodLabels(p) & zoneTable.ClimateZone=="GLOBAL";
    if sum(zoneTable.NCells(zoneTable.Period==periodLabels(p) & ...
              zoneTable.ClimateZone~="GLOBAL")) ~= zoneTable.NCells(v)
        error('Regional counts do not sum to global total in %s.',periodLabels(p));
    end
end

%% ========================= SAVE CSV, XLSX, MAT =========================
writetable(zoneTable, fullfile(outDir,'ClimateZone_GainSummary.csv'));
writetable(winnerTable, fullfile(outDir,'ClimateZone_WinnerShares.csv'));
writetable(categoryTable, fullfile(outDir,'ClimateZone_CategoryLong.csv'));
writetable(validationRecords, fullfile(outDir,'ValidationChecks.csv'));

if writeExcel
    xlsxPath = fullfile(outDir,'ClimateZone_Attribution_Tables.xlsx');
    try
        if isfile(xlsxPath), delete(xlsxPath); end
        writetable(zoneTable,xlsxPath,'Sheet','GainByZone');
        writetable(winnerTable,xlsxPath,'Sheet','WinnerByZone');
        writetable(categoryTable,xlsxPath,'Sheet','CategoryLong');
        writetable(validationRecords,xlsxPath,'Sheet','Validation');
    catch ME
        warning('XLSX export failed (%s). The CSV outputs are complete.',ME.message);
    end
end

save(fullfile(outDir,'ClimateZone_Attribution_Results.mat'), ...
    'lat','lon','validMask','kgFineMap','kgMajorMap', ...
    'climateLetters','climateNames','climateCounts','modelNames', ...
    'years','periodLabels','results','zoneTable','winnerTable', ...
    'categoryTable','validationRecords','matFile','tifFile','-v7.3');

%% ========================= FIGURES (PRIMARY) ===========================
if saveFigures
    zoneRows = zoneTable(zoneTable.Period=="ThreeYearMean" & ...
        zoneTable.ClimateZone~="GLOBAL",:);
    % Rows are ordered A-E.
    if ~isequal(zoneRows.ClimateZone,climateLetters')
        error('Unexpected climate-zone ordering in plot inputs.');
    end

    M = [zoneRows.NoGainPct, zoneRows.LinearPct, zoneRows.NonlinearPct, ...
         zoneRows.HydroclimaticPct, zoneRows.MixedTotalPct];
    f1 = figure('Color','w','Position',[80 80 1000 500]);
    b=bar(1:5,M,'stacked','BarWidth',0.72);
    cmap=[.75 .75 .75; .85 .10 .10; 0 .35 .90; 0 .60 .20; 1 .80 0];
    for j=1:5, b(j).FaceColor=cmap(j,:); end
    ylabel('Share of grid cells within zone (%)');
    xlabel('Koppen-Geiger major climate zone');
    set(gca,'XTick',1:5,'XTickLabel',cellstr(climateLetters),...
        'FontName','Arial','FontSize',12);
    ylim([0 100]); grid on; box on;
    legend({'No gain','Linear','Nonlinear','Hydroclimatic','Mixed'}, ...
        'Location','eastoutside');
    title('Winner-conditioned gain categories by climate zone (2022-2024)');
    print(f1,fullfile(outDir,'Fig1_GainCategories_ByClimate.png'),'-dpng','-r300');
    close(f1);

    W = zeros(5,9);
    for z=1:5
        for m=1:9
            q = winnerTable.Period=="ThreeYearMean" & ...
                winnerTable.ClimateZone==climateLetters(z) & ...
                winnerTable.WinningModel==modelNames(m);
            W(z,m)=winnerTable.PctWithinZone(q);
        end
    end
    f2=figure('Color','w','Position',[80 80 1080 500]);
    bar(1:5,W,'grouped','BarWidth',0.88);
    ylabel('Winner frequency within zone (%)');
    xlabel('Koppen-Geiger major climate zone');
    set(gca,'XTick',1:5,'XTickLabel',cellstr(climateLetters),...
        'FontName','Arial','FontSize',12);
    grid on; box on;
    legend(cellstr(modelNames),'Location','eastoutside');
    title('Best-performing model by climate zone (2022-2024)');
    print(f2,fullfile(outDir,'Fig2_WinningModels_ByClimate.png'),'-dpng','-r300');
    close(f2);

    f3=figure('Color','w','Position',[80 80 850 470]);
    b=bar(1:5,[zoneRows.HistoricalContributionPct, ...
               zoneRows.HydroclimaticContributionPct], 'stacked', 'BarWidth',.65);
    b(1).FaceColor=[.28 .44 .72];
    b(2).FaceColor=[.18 .62 .36];
    ylabel('Share of NRMSE reduction within hydroclimatic-gain cells (%)');
    xlabel('Koppen-Geiger major climate zone');
    set(gca,'XTick',1:5,'XTickLabel',cellstr(climateLetters),...
        'FontName','Arial','FontSize',12);
    ylim([0 100]); grid on; box on;
    legend({'Linear + nonlinear components','Paired hydroclimatic component'}, ...
        'Location','southoutside','Orientation','horizontal');
    title('Within-category gain contributions by climate zone (2022-2024)');
    print(f3,fullfile(outDir,'Fig3_HydroContributions_ByClimate.png'),'-dpng','-r300');
    close(f3);

    f4=figure('Color','w','Position',[80 80 1020 490]);
    h=imagesc(lon, lat, kgMajorMap);
    set(h,'AlphaData',double(validMask));
    set(gca,'YDir','normal','Color',[1 1 1],...
        'FontName','Arial','FontSize',11);
    axis tight; xlabel('Longitude (degrees)');ylabel('Latitude (degrees)');
    colormap([.10 .55 .25; .92 .69 .26; .55 .75 .94; .65 .47 .84; .80 .86 .90]);
    caxis([.5 5.5]);
    cb=colorbar; cb.Ticks=1:5; cb.TickLabels=cellstr(climateLetters);
    title('Koppen-Geiger climate zone assigned to valid GLDAS grid cells');
    print(f4,fullfile(outDir,'Fig4_ClimateZone_GLDASExtent.png'),'-dpng','-r300');
    close(f4);
end

%% ========================= PRINT MAIN SUMMARY ==========================
fprintf('\n========= THREE-YEAR PRIMARY ATTRIBUTION BY CLIMATE =========\n');
disp(zoneTable(zoneTable.Period=="ThreeYearMean", ...
    {'ClimateZone','NCells','LinearPct','NonlinearPct', ...
     'HydroclimaticPct','MixedTotalPct','HistoricalContributionPct', ...
     'HydroclimaticContributionPct'}));
fprintf('Output files saved to:\n%s\n',outDir);

%% ========================= LOCAL FUNCTIONS =============================
function filePath=firstExisting(candidates)
    filePath='';
    for k=1:numel(candidates)
        if isfile(candidates{k}), filePath=candidates{k}; return; end
    end
end

function R=deriveAttribution(E, mask)
% Reconstruct the precise classification / decomposition rules implemented
% by s2_gain_attribution2.m and s3_uc_analysis1.m.
% E: [nLat,nLon,9] NRMSE, model order:
% AR, ARIMA, SARIMA, SVM, RF, LSTM, SVMX, RFX, LSTMX.
    [nLat,nLon,nModel]=size(E);
    if nModel~=9 || ~isequal(size(mask),[nLat,nLon])
        error('Invalid attribution input dimensions.');
    end
    E = double(E);
    for m=1:9
        temp=E(:,:,m);temp(~mask)=NaN;E(:,:,m)=temp;
    end

    [L,li]=min(E(:,:,1:3),[],3);
    [N,ni]=min(E(:,:,4:6),[],3);
    [X,xi]=min(E(:,:,7:9),[],3);
    gStack=cat(3,L,N,X);
    [sE,sIdx]=sort(gStack,3,'ascend');
    bestGroup=sIdx(:,:,1);
    secondGroup=sIdx(:,:,2);
    bestVal=sE(:,:,1);

    winner=nan(nLat,nLon);
    q=mask & bestGroup==1;winner(q)=li(q);
    q=mask & bestGroup==2;winner(q)=ni(q)+3;
    q=mask & bestGroup==3;winner(q)=xi(q)+6;
    winner(~mask)=NaN;

    m1=mask & bestGroup==1 & li==1;
    m2=mask & bestGroup==1 & li~=1;
    m3=mask & bestGroup==2;
    m4=mask & bestGroup==3 & secondGroup==2 & xi==ni;
    m5=mask & bestGroup==3 & secondGroup==1;
    m6=mask & bestGroup==3 & secondGroup==2 & xi~=ni;
    class6=nan(nLat,nLon);
    masks={m1,m2,m3,m4,m5,m6};
    for k=1:6, class6(masks{k})=k; end
    if nnz(mask & ~isfinite(class6))>0
        error('Gain attribution leaves some valid cells unassigned.');
    end
    class5=class6;
    class5(class5==6)=5;  % Mixed-L (5), Mixed-N (6) -> mixed (5)
    class5(class5==5)=5;
    class5(class6==4)=4;

    AR=E(:,:,1);B=bestVal;
    Gtotal=AR-B;
    Gl=zeros(nLat,nLon);Gn=zeros(nLat,nLon);
    Gh=zeros(nLat,nLon);Gm=zeros(nLat,nLon);
    Gl(m2)=AR(m2)-B(m2);
    Gl(m3)=AR(m3)-L(m3); Gn(m3)=L(m3)-N(m3);
    Gl(m4)=AR(m4)-L(m4); Gn(m4)=L(m4)-N(m4);
    Gh(m4)=N(m4)-X(m4);
    Gl(m5)=AR(m5)-L(m5); Gm(m5)=L(m5)-X(m5);
    Gl(m6)=AR(m6)-L(m6); Gn(m6)=L(m6)-N(m6);
    Gm(m6)=N(m6)-X(m6);

    tol=1e-10;
    Gl(abs(Gl)<tol)=0;Gn(abs(Gn)<tol)=0;
    Gh(abs(Gh)<tol)=0;Gm(abs(Gm)<tol)=0;
    Gtotal(abs(Gtotal)<tol)=0;
    reconstructionError=abs(Gl+Gn+Gh+Gm-Gtotal);
    maxErr=max(reconstructionError(mask));
    if maxErr>1e-8
        error('Gain decomposition fails reconstruction check (%.3g).',maxErr);
    end
    Gtotal(~mask)=NaN;
    Gl(~mask)=NaN;Gn(~mask)=NaN;Gh(~mask)=NaN;Gm(~mask)=NaN;

    R=struct('nValid',nnz(mask), ...
        'bestAll9Idx',winner,'gainClassMap',class5,'class6Map',class6, ...
        'AR_NRMSE',AR,'allFusion',B,...
        'G_total',Gtotal,'G_linear',Gl,'G_nonlinear',Gn, ...
        'G_hydro',Gh,'G_mixed',Gm, ...
        'maxReconstructionError',maxErr);
end
