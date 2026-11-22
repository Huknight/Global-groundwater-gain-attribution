function SKLGB_mean_prepare(cfg)
% Construct a compact model master with EXACTLY the same 2000 original cells.
% The original 17,284-cell CommonValid mask and all Noah predictors are
% retained unchanged. Only dGWS_valid changes to the monthly-mean target.
if ~isfile(cfg.originalMasterFile)
    error('Original full-domain GLDAS master is missing: %s',cfg.originalMasterFile);
end
if ~isfile(cfg.spatialTargetFile)
    error('Run SKLGB_mean_run(''upscale'') first: %s',cfg.spatialTargetFile);
end
if ~isfile(cfg.sampleFile)
    error('Fixed sample .mat is missing: %s',cfg.sampleFile);
end
P=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex', ...
    'selectedLat','selectedLon','designWeight','annualNRMSE9', ...
    'annualWinnerModel','annualGainClass','meanGainClass', ...
    'sampleN','populationN','stratumId');
rows=double(P.selectedValidRow(:));ix=double(P.selectedLinearIndex(:));
assert(P.sampleN==2000 && P.populationN==17284);
assert(numel(rows)==2000 && numel(unique(rows))==2000 && ...
    numel(unique(ix))==2000 && all(diff(rows)>0) && ...
    all(diff(ix)>0) && abs(sum(double(P.designWeight(:)))-17284)<1e-6, ...
    '2000-grid sample or strata weights invalid.');
if isfile(cfg.pilotMasterFile) && ~cfg.rebuildPilotMaster
    S=load(cfg.pilotMasterFile,'Meta','originalValidRows','validLinearIndex');
    if ~isfield(S,'Meta') || ...
            ~strcmp(S.Meta.targetDefinition,'difference of consecutive 1-degree monthly-mean GWS states') || ...
            ~isfield(S.Meta,'requireFullDailyCoverage') || ...
            S.Meta.requireFullDailyCoverage~=cfg.requireFullDailyCoverage || ...
            ~isequal(double(S.originalValidRows(:)),rows) || ...
            ~isequal(double(S.validLinearIndex(:)),ix)
        error('Previously prepared sample is incompatible. Inspect before setting rebuildPilotMaster=true.');
    end
    fprintf('Reusing existing monthly-mean 2000 master: %s\n',cfg.pilotMasterFile);
    return;
end
D=load(cfg.originalMasterFile,'dGWS_valid', 'Prec_valid','Evap_valid', ...
    'Runoff_valid','dSM2m_valid','SurfT_valid','CommonMask', ...
    'CommonValid','validLinearIndex','validLat','validLon','lat','lon','ym');
assert(numel(D.validLinearIndex)==17284 && ...
    isequal(size(D.CommonValid),[150,360]) && ...
    isequal(double(D.validLinearIndex(rows(:))),ix) && ...
    isequal(double(D.validLat(rows)),double(P.selectedLat(:))) && ...
    isequal(double(D.validLon(rows)),double(P.selectedLon(:))), ...
    'Original sample identities/coordinates differ from the reference GLDAS mask.');
R=load(cfg.spatialTargetFile,'dGWS_mm','lat','lon','ym','Meta');
assert(isequal(size(R.dGWS_mm),[150,360,252]),'Alternative target must be 150x360x252.');
assert(isequal(double(R.ym(:)),double(D.ym(:))), ...
    'Alternative dGWS months do not align with the original time axis.');
assert(all(abs(double(R.lat(:))-double(D.lat(:)))<1e-8) && ...
    all(abs(double(R.lon(:))-double(D.lon(:)))<1e-8), ...
    'Alternative and original 1-degree grids are not aligned.');
assert(strcmp(R.Meta.targetDefinition, ...
    'difference of consecutive 1-degree monthly-mean GWS states'));
Yall=reshape(double(R.dGWS_mm),150*360,252);
dGWS_valid=Yall(ix,:);
if any(~isfinite(dGWS_valid(:)))
    nBad=nnz(any(~isfinite(dGWS_valid),2));
    error(['Alternative target has missing/invalid months in %d of 2000 ' ...
           'original selected cells. No replacement cells or silent masking allowed.'],nBad);
end
originalMaskFull=logical(D.CommonValid);
nMaskComplete=sum(all(isfinite(Yall(double(D.validLinearIndex(:)),:)),2));
if nMaskComplete~=17284
    warning('SKLGB:MonthlyMeanCoverage', ...
        '%d/%d original global mask cells complete under target; 2000 selected cells are checked separately.', ...
        nMaskComplete,17284);
end
Prec_valid=double(D.Prec_valid(rows,:));
Evap_valid=double(D.Evap_valid(rows,:));
Runoff_valid=double(D.Runoff_valid(rows,:));
dSM2m_valid=double(D.dSM2m_valid(rows,:));
SurfT_valid=double(D.SurfT_valid(rows,:));
lat=double(D.lat(:));lon=double(D.lon(:));ym=double(D.ym(:));
validLinearIndex=ix(:);
validLat=double(P.selectedLat(:));validLon=double(P.selectedLon(:));
CommonValid=false(150,360);CommonValid(ix)=true;
CommonMask=nan(150,360);CommonMask(ix)=1;
originalValidRows=rows(:);
pilotDesignWeight=double(P.designWeight(:));
pilotStratumId=double(P.stratumId(:));
pilotBaselineNRMSE9=double(P.annualNRMSE9);
pilotMeanGainClass=double(P.meanGainClass(:));
pilotAnnualGainClass=double(P.annualGainClass);
pilotAnnualWinnerModel=double(P.annualWinnerModel);
selectedSamplePositions=(1:2000)';
inputVars={dGWS_valid,Prec_valid,Evap_valid,Runoff_valid,dSM2m_valid,SurfT_valid};
for k=1:6
    assert(isequal(size(inputVars{k}),[2000,252]) && ...
        all(isfinite(inputVars{k}(:))), 'Predictor %d incomplete.',k);
end
assert(isequal(find(CommonValid),ix(:)));
[Lon2D,Lat2D]=meshgrid(lon,lat);
assert(all(abs(Lat2D(ix)-validLat)<1e-8) && ...
    all(abs(Lon2D(ix)-validLon)<1e-8));
% Prove GLDAS-Noah predictors are copied, not recalculated.
Meta=struct();
Meta.targetDefinition=R.Meta.targetDefinition;
Meta.sampleN=2000;Meta.populationN=17284;
Meta.originalMaskFullCells=nnz(originalMaskFull);
Meta.completeOriginalMaskCellsUnderAlternativeTarget=nMaskComplete;
Meta.noahVariables='Unchanged: Prec, Evap, Runoff, dSM2m, SurfT';
Meta.originalFullMaster=cfg.originalMasterFile;
Meta.targetSource=cfg.spatialTargetFile;
Meta.testYears=[2022 2023 2024];
Meta.requireFullDailyCoverage=R.Meta.requireFullDailyCoverage;
Meta.noPopulationExtrapolation=true;
if ~isfolder(fileparts(cfg.pilotMasterFile)),mkdir(fileparts(cfg.pilotMasterFile));end
f=[cfg.pilotMasterFile,'.tmp.mat'];
save(f,'dGWS_valid','Prec_valid','Evap_valid','Runoff_valid', ...
    'dSM2m_valid','SurfT_valid','CommonValid','CommonMask', ...
    'validLinearIndex','validLat','validLon','lat','lon','ym', ...
    'originalValidRows','selectedSamplePositions','pilotBaselineNRMSE9', ...
    'pilotMeanGainClass','pilotAnnualGainClass','pilotAnnualWinnerModel', ...
    'pilotDesignWeight','pilotStratumId','Meta','-v7.3');
[ok,msg]=movefile(f,cfg.pilotMasterFile,'f');
if ~ok,error('Failed to save selected2000 master: %s',msg);end
fprintf('Prepared 2000 identical grid cells; same five Noah predictors.\n');
fprintf('Original full mask complete under new target: %d / 17284.\n',nMaskComplete);
end
