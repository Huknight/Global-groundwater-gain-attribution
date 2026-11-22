function SKLGB_ext_run(action)
% Call SKLGB_ext_run('preflight') first, then SKLGB_ext_run('train').
% Call SKLGB_ext_run('analyze') after all nine models complete.
if nargin<1, action='train'; end
if ~ismember(action,{'preflight','train','analyze'})
    error('action must be preflight, train, or analyze.');
end
cfg=SKLGB_ext_config();
here=fileparts(mfilename('fullpath'));
addpath(here);
if ~isfolder(cfg.outputRoot), mkdir(cfg.outputRoot); end
if strcmp(action,'analyze')
    SKLGB_ext_analyze(cfg);
    return;
end
if ~isfile(cfg.sampleFile), error('Missing fixed 2000 sample MAT file.'); end
if ~isfile(cfg.originalMasterFile)
    error('Missing original GLDAS master: %s',cfg.originalMasterFile);
end
P=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex', ...
    'designWeight','annualNRMSE9','sampleN','populationN');
assert(numel(P.selectedValidRow)==2000 && ...
    numel(unique(P.selectedLinearIndex))==2000 && ...
    P.sampleN==2000 && P.populationN==17284, ...
    'Sample identity/population-size mismatch.');
assert(all(isfinite(P.annualNRMSE9(:))) && ...
    abs(sum(double(P.designWeight))-17284)<1e-6, ...
    'Sample reference NRMSE or stratification weights invalid.');
SKLGB_ext_prepare_master(cfg);
D=load(cfg.pilotMasterFile,'originalValidRows','originalLinearIndex', ...
    'pilotDesignWeight','dGWS_valid','ym');
assert(isequal(double(P.selectedValidRow(:)),double(D.originalValidRows(:))));
assert(isequal(double(P.selectedLinearIndex(:)),double(D.originalLinearIndex(:))));
assert(size(D.dGWS_valid,1)==2000 && size(D.dGWS_valid,2)==252);
assert(numel(D.ym)==252 && D.ym(1)==200401 && D.ym(end)==202412);
fprintf('\n=======================================================\n');
fprintf('SKLGB: extended test years 2019, 2020, 2021\n');
fprintf('SINGLE PC, fixed 2000 cells / original 17,284-row IDs\n');
fprintf('Original RNG seed: 20260605; original hyperparameters\n');
fprintf('Per T: train 2004..T-2, validate T-1, test T\n');
fprintf('Outputs: %s\n',cfg.outputRoot);
fprintf('=======================================================\n');
if strcmp(action,'preflight')
    fprintf('PREFLIGHT PASSED; models scheduled: %s\n',strjoin(cfg.modelsToRun,', '));
    return;
end
if cfg.useParallel
    pool=gcp('nocreate');
    if isempty(pool), parpool('local',cfg.numWorkers); end
end
models=cfg.modelsToRun;
for j=1:numel(models)
    mdl=models{j};
    name=fullfile(cfg.outputRoot,mdl,sprintf('%s_extended_2019_2021.mat',mdl));
    if isfile(name)
        A=load(name,'Meta','testYears','validLinearIndex');
        if ~isequal(double(A.validLinearIndex(:)), ...
                double(D.originalLinearIndex(:))) || ...
                ~isequal(double(A.testYears(:))',[2019 2020 2021]) || ...
                A.Meta.randomSeed~=20260605
            error('Existing output is incompatible: %s',name);
        end
        if ismember(mdl,{'AR','ARIMA','SARIMA'})
            if ~isfield(A.Meta,'provenance') || ...
                ~strcmp(A.Meta.provenance,'ORIGINAL_USER_SUPPLIED_MATLAB_SOURCE') || ...
                ~isfield(A.Meta,'sourceCoreSHA256')
                error(['Existing %s result is from the old RECONSTRUCTED ' ...
                    'linear implementation. It will NOT be silently reused. ' ...
                    'Back up/remove only that result file, then restart: %s'],mdl,name);
            end
            src=fullfile(here,['SKLGB_ext_' mdl '_core.m']);
            md=java.security.MessageDigest.getInstance('SHA-256');
            fid=fopen(src,'rb');
            if fid<0,error('Missing source core: %s',src);end
            cc=onCleanup(@() fclose(fid));
            while ~feof(fid)
                chunk=fread(fid,65536,'*uint8');
                if ~isempty(chunk),md.update(typecast(chunk,'int8'));end
            end
            b=typecast(int8(md.digest()),'uint8');
            existingHash=lower(reshape(dec2hex(b,2).',1,[]));
            clear cc
            if ~strcmp(existingHash,A.Meta.sourceCoreSHA256)
                error('Previously completed %s output used a different source core: %s', ...
                    mdl,name);
            end
        end
        fprintf('[%d/%d] %s verified existing -> SKIP\n',j,numel(models),mdl);
        continue;
    end
    fprintf('\n[%d/%d] TRAIN %s\n',j,numel(models),mdl);
    ticModel=tic;
    if ismember(mdl,{'AR','ARIMA','SARIMA'})
        SKLGB_ext_linear(cfg,mdl);
    else
        feval(['SKLGB_ext_' mdl],cfg);
    end
    A=load(name,'Meta','testYears','validLinearIndex', ...
        [mdl '_NRMSE_valid'],[mdl '_TestSuccess_valid']);
    ok=A.([mdl '_TestSuccess_valid']);
    fprintf('%s: success %d/%d grid-year forecasts, wall %.1f min.\n', ...
        mdl,nnz(ok),numel(ok),toc(ticModel)/60);
    logFile=fullfile(cfg.outputRoot,'extended_run_log.csv');
    first=~isfile(logFile);
    fid=fopen(logFile,'a');
    if fid<0,error('Cannot open execution log');end
    if first, fprintf(fid,'Model,Status,SuccessGridYears,TotalGridYears,Minutes\n');end
    fprintf(fid,'%s,completed,%d,%d,%.6f\n', ...
        mdl,nnz(ok),numel(ok),toc(ticModel)/60);
    fclose(fid);
end
fprintf('\nAll scheduled models finished.\n');
required={'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};
allHere=true;
for k=1:numel(required)
    f=fullfile(cfg.outputRoot,required{k}, ...
        sprintf('%s_extended_2019_2021.mat',required{k}));
    allHere=allHere && isfile(f);
end
if allHere
    fprintf('All nine model files exist; running annual summaries.\n');
    SKLGB_ext_analyze(cfg);
else
    fprintf('Some model outputs are pending; call SKLGB_ext_run(''analyze'') after all nine finish.\n');
end
end
