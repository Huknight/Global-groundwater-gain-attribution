function SKLGB_mean_extract(cfg, mode)
% Read ALL 7702 daily GWS_tavg records, 2003-12-01--2024-12-31.
% Resumable by month. DO NOT modify the original GLDAS NetCDF files.
if nargin<2, mode='extract'; end
if ~ismember(mode,{'preflight','extract'})
    error('Mode must be preflight or extract.');
end
if ~isfolder(cfg.rawDAFolder)
    error('Raw DA daily folder not found: %s',cfg.rawDAFolder);
end
months=(datetime(2003,12,1):calmonths(1):datetime(2024,12,1))';
nMonth=numel(months);
assert(nMonth==253);
firstFile=dailyFile(cfg,datetime(2003,12,1));
if ~isfile(firstFile),error('Missing first daily GWS source: %s',firstFile);end
latRaw=double(ncread(firstFile,'lat'));
lonRaw=double(ncread(firstFile,'lon'));
lonRaw(lonRaw>180)=lonRaw(lonRaw>180)-360;
[lat,latOrder]=sort(latRaw(:),'ascend');
[lon,lonOrder]=sort(lonRaw(:),'ascend');
assert(numel(lat)==600 && numel(lon)==1440 && ...
    all(abs(diff(lat)-0.25)<1e-5) && all(abs(diff(lon)-0.25)<1e-5), ...
    'Unexpected GLDAS-DA native grid.');
assert(abs(lat(1)+59.875)<1e-5 && abs(lat(end)-89.875)<1e-5 && ...
    abs(lon(1)+179.875)<1e-5 && abs(lon(end)-179.875)<1e-5, ...
    'Unexpected native grid orientation/extent.');
info=ncinfo(firstFile);
assert(any(strcmp({info.Variables.Name},'GWS_tavg')),'GWS_tavg variable is missing.');
if strcmp(mode,'preflight')
    if cfg.checkRawFilesAtPreflight
        missing=cell(0,1); checked=0;
        for j=1:nMonth
            ds=months(j)+days((0:eomday(year(months(j)),month(months(j)))-1)');
            for k=1:numel(ds)
                f=dailyFile(cfg,ds(k));
                if ~isfile(f)
                    missing{end+1,1}=f; %#ok<AGROW>
                end
                checked=checked+1;
            end
        end
        if ~isempty(missing)
            fprintf('First missing files:\n');
            for k=1:min(10,numel(missing)),fprintf('%s\n',missing{k});end
            error('%d of %d mandatory daily GLDAS files are missing.',numel(missing),checked);
        end
        assert(checked==7702,'Expected 7702 daily files; got %d.',checked);
        fprintf('PREFLIGHT: all %d daily files found, 253 months.\n',checked);
    end
    fprintf('DA first file validated: 600x1440, GWS_tavg, native coordinates OK.\n');
    return;
end
if ~isfolder(cfg.nativeMonthDir),mkdir(cfg.nativeMonthDir);end
fprintf('Extracting %d monthly-mean states (%d raw daily files).\n',nMonth,7702);
if cfg.useParallel
    pool=gcp('nocreate');
    if isempty(pool),parpool('local',cfg.numWorkers);
    elseif pool.NumWorkers~=cfg.numWorkers
        error('Pool has %d workers but cfg.numWorkers=%d.',pool.NumWorkers,cfg.numWorkers);
    end
    % Parfor processes separate months with independent output files.
    parfor j=1:nMonth
        processMonth(cfg,months(j),latRaw,lonRaw,latOrder,lonOrder,lat,lon);
    end
else
    for j=1:nMonth
        fprintf('Month %d/%d: %s\n',j,nMonth,char(months(j)));
        processMonth(cfg,months(j),latRaw,lonRaw,latOrder,lonOrder,lat,lon);
    end
end
fprintf('All 253 monthly-mean native states have been computed or verified.\n');
end

function processMonth(cfg,m,latRaw,lonRaw,latOrder,lonOrder,lat,lon)
monthID=year(m)*100+month(m);
nDays=eomday(year(m),month(m));
outFile=fullfile(cfg.nativeMonthDir,sprintf('GWS_native_monthmean_%06d.mat',monthID));
if isfile(outFile)
    Q=load(outFile,'Meta');
    if ~isfield(Q,'Meta') || ~isfield(Q.Meta,'monthID') || ...
            Q.Meta.monthID~=monthID || Q.Meta.days~=nDays || ...
            Q.Meta.requireFullDailyCoverage~=cfg.requireFullDailyCoverage || ...
            ~isfield(Q.Meta,'procedureVersion') || Q.Meta.procedureVersion~=1
        error('Existing monthly checkpoint incompatible: %s',outFile);
    end
    return;
end
nLat=numel(lat);nLon=numel(lon);
S=zeros(nLat,nLon); counts=zeros(nLat,nLon,'uint8');
for k=1:nDays
    d=datetime(year(m),month(m),k);
    nc=dailyFile(cfg,d);
    if ~isfile(nc),error('Missing daily GLDAS-DA source: %s',nc);end
    if k==1 || k==nDays || cfg.checkEveryDailyCoordinates
        checkCoordinates(nc,latRaw,lonRaw);
    end
    A=readGWS(nc,nLat,nLon,latOrder,lonOrder);
    ok=isfinite(A);
    S(ok)=S(ok)+A(ok);
    counts(ok)=counts(ok)+uint8(1);
end
GWS_mean025=S./double(counts);
GWS_mean025(counts==0)=NaN;
if cfg.requireFullDailyCoverage
    GWS_mean025(counts<nDays)=NaN;
end
Meta=struct('monthID',monthID,'days',nDays, ...
    'validFullDays',nnz(counts==nDays),'partialDays',nnz(counts>0 & counts<nDays), ...
    'noDays',nnz(counts==0),'requireFullDailyCoverage', ...
    cfg.requireFullDailyCoverage,'procedureVersion',1, ...
    'source','GLDAS_CLSM025_DA1_D.AYYYYMMDD.022.nc4', ...
    'variable','GWS_tavg','units','mm','latFirst',lat(1),'lonFirst',lon(1));
tmp=[outFile,'.tmp.mat'];
save(tmp,'GWS_mean025','Meta','-v7.3');
[ok,msg]=movefile(tmp,outFile,'f');
if ~ok,error('Failed to save month %d: %s',monthID,msg);end
end

function nc=dailyFile(cfg,d)
y=year(d); md=year(d)*10000+month(d)*100+day(d);
base=sprintf('GLDAS_CLSM025_DA1_D.A%08d.022',md);
ydir=fullfile(cfg.rawDAFolder,sprintf('%04d',y));
nc=fullfile(ydir,[base '.nc4']);
if ~isfile(nc)
    alt=fullfile(ydir,base);
    if isfile(alt),nc=alt;end
end
end

function checkCoordinates(nc,latRaw,lonRaw)
a=double(ncread(nc,'lat'));b=double(ncread(nc,'lon'));
b(b>180)=b(b>180)-360;
if numel(a)~=numel(latRaw) || numel(b)~=numel(lonRaw) || ...
        any(abs(a(:)-latRaw(:))>1e-5) || any(abs(b(:)-lonRaw(:))>1e-5)
    error('Daily native coordinate mismatch: %s',nc);
end
end

function A=readGWS(nc,nLat,nLon,latOrder,lonOrder)
raw=squeeze(double(ncread(nc,'GWS_tavg')));
fills=[];
try, x=double(ncreadatt(nc,'GWS_tavg','_FillValue'));fills=[fills;x(:)];catch,end
try, x=double(ncreadatt(nc,'GWS_tavg','missing_value'));fills=[fills;x(:)];catch,end
for k=1:numel(fills),raw(raw==fills(k))=NaN;end
raw(~isfinite(raw) | abs(raw)>1e19)=NaN;
if isequal(size(raw),[nLon,nLat]),A=raw';
elseif isequal(size(raw),[nLat,nLon]),A=raw;
else,error('Unexpected GWS_tavg size in %s: %s',nc,mat2str(size(raw)));
end
A=A(latOrder,lonOrder);
end
