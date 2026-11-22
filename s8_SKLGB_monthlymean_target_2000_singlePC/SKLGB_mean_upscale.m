function SKLGB_mean_upscale(cfg)
% Native monthly means -> spherical-area-weighted 1-degree monthly means.
% Target EXACTLY matches the manuscript's monthly-mean difference definition:
%   dGWS1(:,:,t) = GWS1_mean(:,:,t+1) - GWS1_mean(:,:,t).
% A separate native-paired aggregation is computed for diagnostic auditing.
months=(datetime(2003,12,1):calmonths(1):datetime(2024,12,1))';
assert(numel(months)==253);
if isfile(cfg.spatialTargetFile)
    Q=load(cfg.spatialTargetFile,'Meta');
    if isfield(Q,'Meta') && isfield(Q.Meta,'targetDefinition') && ...
            strcmp(Q.Meta.targetDefinition,'difference of consecutive 1-degree monthly-mean GWS states') && ...
            isfield(Q.Meta,'requireFullDailyCoverage') && ...
            Q.Meta.requireFullDailyCoverage==cfg.requireFullDailyCoverage
        fprintf('Existing upscaled monthly-mean target verified: %s\n',cfg.spatialTargetFile);
        return;
    end
    error('Alternative target already exists but does not match this procedure: %s',cfg.spatialTargetFile);
end
f0=monthlyFile(cfg,months(1));
if ~isfile(f0),error('Missing December 2003 monthly mean: %s',f0);end
GWS1_mm=nan(150,360,253);
ValidAreaFraction=struct();
ValidAreaFraction.GWS1_mm=nan(150,360,253,'single');
ValidAreaFraction.dGWS_mm=nan(150,360,252,'single');
ValidAreaFraction.nativePaired_mm=nan(150,360,252,'single');
dGWS_nativePaired_mm=nan(150,360,252);
lat025=(-59.875:0.25:89.875)';
lon025=(-179.875:0.25:179.875)';
lat=mean(reshape(lat025,4,150),1)';
lon=mean(reshape(lon025,4,360),1)';
[Lon2D,Lat2D]=meshgrid(lon,lat);
w=sind(lat025+0.125)-sind(lat025-0.125);
assert(all(w>0));
prev=[]; nPartial=0; nEmpty=0;
for j=1:numel(months)
    f=monthlyFile(cfg,months(j));
    if ~isfile(f),error('Missing month %s. Run extract first.',char(months(j)));end
    Q=load(f,'GWS_mean025','Meta');
    assert(Q.Meta.monthID==year(months(j))*100+month(months(j)));
    assert(Q.Meta.requireFullDailyCoverage==cfg.requireFullDailyCoverage);
    A=double(Q.GWS_mean025);
    assert(isequal(size(A),[600,1440]));
    [g,c]=weighted025to1(A,w);
    GWS1_mm(:,:,j)=g;
    ValidAreaFraction.GWS1_mm(:,:,j)=c;
    nPartial=nPartial+Q.Meta.partialDays;
    nEmpty=nEmpty+Q.Meta.noDays;
    if j>1
        % Paired-cell diagnostic follows the original Step 4 spatial rule.
        [paired,pc]=weighted025to1(A-prev,w);
        dGWS_nativePaired_mm(:,:,j-1)=paired;
        ValidAreaFraction.nativePaired_mm(:,:,j-1)=pc;
    end
    prev=A;
    if mod(j,12)==0 || j==253
        fprintf('Aggregated %3d/253 monthly GWS states.\n',j);
    end
end
% Article definition: the difference of consecutive monthly MEAN 1-degree GWS.
dGWS_mm=diff(GWS1_mm,1,3);
ValidAreaFraction.dGWS_mm=min(ValidAreaFraction.GWS1_mm(:,:,1:252), ...
    ValidAreaFraction.GWS1_mm(:,:,2:253));
ymState=double(year(months)*100+month(months));
ym=ymState(2:end);yearVec=floor(ym/100);monthVec=mod(ym,100);
finiteBoth=isfinite(dGWS_mm)&isfinite(dGWS_nativePaired_mm);
if any(finiteBoth(:))
    dd=abs(dGWS_mm(finiteBoth)-dGWS_nativePaired_mm(finiteBoth));
    CompareAudit=struct('nComparable',numel(dd),'meanAbsDiff',mean(dd), ...
        'maxAbsDiff',max(dd),'nDiffOver1e8',nnz(dd>1e-8), ...
        'nDiffOver1e4',nnz(dd>1e-4), ...
        'comment','Source-cell coverage variations cause differences between native-paired change and change of 1-degree monthly means.');
else
    CompareAudit=struct('nComparable',0,'meanAbsDiff',NaN,'maxAbsDiff',NaN,...
        'nDiffOver1e8',0,'nDiffOver1e4',0);
end
Meta=struct();
Meta.sourceProduct='GLDAS_CLSM025_DA1_D_2.2';
Meta.variable='GWS_tavg';Meta.sourceDailyPeriod='2003-12-01 to 2024-12-31';
Meta.statePeriod='2003-12 to 2024-12';Meta.targetPeriod='2004-01 to 2024-12';
Meta.targetDefinition='difference of consecutive 1-degree monthly-mean GWS states';
Meta.spatialUpscaling='Same spherical-area weights as original Step 4, separately per monthly-mean state';
Meta.weightExpression='sind(latitude+0.125)-sind(latitude-0.125)';
Meta.requireFullDailyCoverage=cfg.requireFullDailyCoverage;
Meta.nativePixelPartialDayCases=nPartial;
Meta.nativePixelNoDayCases=nEmpty;
Meta.note='Monthly means use daily GWS_tavg, NOT the last day of each month.';
Meta.compareNativePaired=CompareAudit;
Meta.nDays=7702;
if ~isfolder(fileparts(cfg.spatialTargetFile))
    mkdir(fileparts(cfg.spatialTargetFile));
end
fprintf('Monthly GWS: %s\n',mat2str(size(GWS1_mm)));
fprintf('Alternative target: %s; paired diagnostic max difference: %.6g mm\n', ...
    mat2str(size(dGWS_mm)),CompareAudit.maxAbsDiff);
tmp=[cfg.spatialTargetFile,'.tmp.mat'];
save(tmp,'GWS1_mm','dGWS_mm','dGWS_nativePaired_mm', ...
    'ValidAreaFraction','lat','lon','Lat2D','Lon2D','ym','ymState', ...
    'yearVec','monthVec','Meta','CompareAudit','-v7.3');
[ok,msg]=movefile(tmp,cfg.spatialTargetFile,'f');
if ~ok,error('Failed final spatial-target save: %s',msg);end
end
function f=monthlyFile(cfg,m)
    f=fullfile(cfg.nativeMonthDir,...
        sprintf('GWS_native_monthmean_%06d.mat',year(m)*100+month(m)));
end
function [Y,C]=weighted025to1(X,wLat)
    assert(isequal(size(X),[600,1440]) && numel(wLat)==600);
    X=reshape(X,4,150,4,360);w=reshape(wLat,4,150,1,1);
    valid=isfinite(X);X(~valid)=0;
    num=reshape(sum(sum(X.*w,1),3),150,360);
    den=reshape(sum(sum(double(valid).*w,1),3),150,360);
    Y=num./den;Y(den==0)=NaN;
    denomFull=4*sum(reshape(wLat,4,150),1)';
    C=single(den./denomFull);
    C=max(0,min(1,C));
end
