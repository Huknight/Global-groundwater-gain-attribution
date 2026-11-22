function SKLGB_mean_smoke(model)
% Run the original model core on 10 selected monthly-mean-target grid cells.
% Uses a SEPARATE smoke directory and cannot alter full-scale checkpoints.
if nargin<1,model='AR';end
allowed={'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};
if ~ismember(model,allowed),error('Unknown model: %s',model);end
cfg=SKLGB_mean_config();
addpath(fileparts(mfilename('fullpath')));
if ~isfile(cfg.pilotMasterFile)
    error('Prepare alternative target and sample master before smoke.');
end
P=load(cfg.sampleFile,'nestedPriority');
M=load(cfg.pilotMasterFile);
rows=find(double(P.nestedPriority(:))<=10);
assert(numel(rows)==10);
vars={'dGWS_valid','Prec_valid','Evap_valid','Runoff_valid', ...
    'dSM2m_valid','SurfT_valid'};
for k=1:numel(vars),M.(vars{k})=M.(vars{k})(rows,:);end
vars={'originalValidRows','validLinearIndex','validLat','validLon', ...
    'pilotDesignWeight','pilotStratumId','pilotMeanGainClass'};
for k=1:numel(vars),M.(vars{k})=M.(vars{k})(rows,:);end
M.pilotBaselineNRMSE9=M.pilotBaselineNRMSE9(rows,:,:);
M.pilotAnnualGainClass=M.pilotAnnualGainClass(rows,:);
M.pilotAnnualWinnerModel=M.pilotAnnualWinnerModel(rows,:);
M.selectedSamplePositions=M.selectedSamplePositions(rows,:);
M.CommonValid=false(150,360);M.CommonValid(M.validLinearIndex)=true;
M.CommonMask=nan(150,360);M.CommonMask(M.validLinearIndex)=1;
cfg.outputRoot=fullfile(cfg.outputRoot,'smoke10');
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
cfg.pilotMasterFile=fullfile(cfg.outputRoot,'smoke10_monthmean_master.mat');
save(cfg.pilotMasterFile,'-struct','M','-v7.3');
cfg.useParallel=false;cfg.numWorkers=1;cfg.block_num=1;
cfg.nLocal=10;
f=fullfile(cfg.outputRoot,model,sprintf('%s_monthmean_2022_2024.mat',model));
if isfile(f),delete(f);end
% Clear smoke-only checkpoint to ensure that each smoke is a real fresh run.
checkpoint=fullfile(cfg.outputRoot,model, ...
    sprintf('%s_monthmean_2022_2024_checkpoint.mat',model));
if isfile(checkpoint),delete(checkpoint);end
linearCheckpoint=fullfile(cfg.outputRoot,model, ...
    sprintf('%s_monthmean_2022_2024_originalcore_v2_checkpoint.mat',model));
if isfile(linearCheckpoint),delete(linearCheckpoint);end
if ismember(model,{'AR','ARIMA','SARIMA'})
    SKLGB_mean_linear(cfg,model);
else
    feval(['SKLGB_mean_' model],cfg);
end
S=load(f,[model '_TestSuccess_valid']);
ok=S.([model '_TestSuccess_valid']);
fprintf('SMOKE %s: %d/%d successful grid-year forecasts.\n', ...
    model,nnz(ok),numel(ok));
if nnz(ok)~=numel(ok)
    warning('SKLGB:SmokeFailure','Some %s smoke forecasts failed; inspect failure codes.',model);
end
end
