function SKLGB_ext_smoke(model)
% Optional small-sample execution test; never overwrites the full 2000 run.
% Example: SKLGB_ext_smoke('RF')
if nargin<1,model='RF';end
if ~ismember(model,{'AR','ARIMA','SARIMA','SVM','RF','LSTM', ...
        'SVMX','RFX','LSTMX'})
    error('Unknown model');
end
cfg=SKLGB_ext_config();
addpath(fileparts(mfilename('fullpath')));
if ~isfile(cfg.pilotMasterFile)
    error('Run SKLGB_ext_run(''preflight'') first.');
end
P=load(cfg.sampleFile,'nestedPriority');
M=load(cfg.pilotMasterFile);
rows=find(P.nestedPriority(:)<=10);
assert(numel(rows)==10);
vars={'dGWS_valid','Prec_valid','Evap_valid','Runoff_valid', ...
    'dSM2m_valid','SurfT_valid'};
for j=1:numel(vars)
    M.(vars{j})=M.(vars{j})(rows,:);
end
vars={'originalValidRows','originalLinearIndex','validLinearIndex', ...
    'validLat','validLon','pilotDesignWeight','pilotStratumId', ...
    'pilotMeanGainClass'};
for j=1:numel(vars)
    M.(vars{j})=M.(vars{j})(rows,:);
end
M.pilotBaselineNRMSE9=M.pilotBaselineNRMSE9(rows,:,:);
M.pilotAnnualGainClass=M.pilotAnnualGainClass(rows,:);
M.pilotAnnualWinnerModel=M.pilotAnnualWinnerModel(rows,:);
M.selectedSamplePositions=M.selectedSamplePositions(rows,:);
M.CommonValid=false(150,360);
M.CommonValid(M.validLinearIndex)=true;
M.CommonMask=nan(150,360);
M.CommonMask(M.validLinearIndex)=1;
cfg.outputRoot=fullfile(cfg.outputRoot,'smoke10');
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
cfg.pilotMasterFile=fullfile(cfg.outputRoot,'smoke10_master.mat');
save(cfg.pilotMasterFile,'-struct','M','-v7.3');
cfg.useParallel=false;
cfg.numWorkers=1;
cfg.block_num=1;
cfg.nLocal=10;
out=fullfile(cfg.outputRoot,model, ...
    sprintf('%s_extended_2019_2021.mat',model));
if isfile(out),delete(out);end
switch model
    case {'AR','ARIMA','SARIMA'}
        SKLGB_ext_linear(cfg,model);
    otherwise
        feval(['SKLGB_ext_' model],cfg);
end
S=load(out,[model '_TestSuccess_valid']);
ok=S.([model '_TestSuccess_valid']);
fprintf('SMOKE %s: %d/%d successful grid-year tests.\n', ...
    model,nnz(ok),numel(ok));
if nnz(ok)<numel(ok)
    warning('Some smoke tests failed; inspect the model file before full run.');
end
end
