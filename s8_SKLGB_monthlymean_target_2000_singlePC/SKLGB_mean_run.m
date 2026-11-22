function SKLGB_mean_run(action)
% Monthly mean GWS target sensitivity; fixed 2000 cells, original nine models.
% Actions: preflight, extract, upscale, prepare, preprocess, train, analyze, all.
if nargin<1,action='preflight';end
allowed={'preflight','extract','upscale','prepare','preprocess','train','analyze','all'};
if ~ismember(action,allowed),error('Unknown action: %s',action);end
cfg=SKLGB_mean_config();
here=fileparts(mfilename('fullpath'));addpath(here);
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
if ismember(action,{'preflight','all'})
    if ~isfile(cfg.originalMasterFile)
        error('Original global GLDAS master missing: %s',cfg.originalMasterFile);
    end
    if ~isfile(cfg.sampleFile)
        error('Original 2000-grid reference file missing: %s',cfg.sampleFile);
    end
    Q=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex', ...
        'sampleN','populationN','designWeight','annualNRMSE9');
    assert(Q.sampleN==2000 && Q.populationN==17284 && ...
        numel(unique(double(Q.selectedLinearIndex)))==2000 && ...
        abs(sum(double(Q.designWeight(:)))-17284)<1e-6 && ...
        isequal(size(Q.annualNRMSE9),[2000,9,3]), ...
        'Fixed sample identity, weighting or reference result mismatch.');
    SKLGB_mean_extract(cfg,'preflight');
    fprintf('PREFLIGHT PASSED. Original global mask and fixed sample recognized.\n');
    if strcmp(action,'preflight'),return;end
end
if ismember(action,{'extract','preprocess','all'})
    SKLGB_mean_extract(cfg,'extract');
    if strcmp(action,'extract'),return;end
end
if ismember(action,{'upscale','preprocess','all'})
    SKLGB_mean_upscale(cfg);
    if strcmp(action,'upscale'),return;end
end
if ismember(action,{'prepare','preprocess','all'})
    SKLGB_mean_prepare(cfg);
    if ismember(action,{'prepare','preprocess'}),return;end
end
if strcmp(action,'analyze')
    SKLGB_mean_analyze(cfg);
    return;
end
if ismember(action,{'train','all'})
    if ~isfile(cfg.pilotMasterFile)
        error('No monthly-mean 2000-cell master. Run SKLGB_mean_run(''preprocess'') first.');
    end
    D=load(cfg.pilotMasterFile,'dGWS_valid','originalValidRows', ...
        'validLinearIndex','Meta');
    P=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex');
    assert(isequal(size(D.dGWS_valid),[2000,252]) && ...
        isequal(double(D.originalValidRows(:)),double(P.selectedValidRow(:))) && ...
        isequal(double(D.validLinearIndex(:)),double(P.selectedLinearIndex(:))) && ...
        strcmp(D.Meta.targetDefinition, ...
         'difference of consecutive 1-degree monthly-mean GWS states'), ...
        'Alternative-target model master not consistent with original sample.');
    fprintf('MONTH-MEAN SENSITIVITY: 2000 fixed cells, 2022-2024, original parameter settings.\n');
    models=cfg.modelsToRun;
    for j=1:numel(models)
        model=models{j};
        f=fullfile(cfg.outputRoot,model,sprintf('%s_monthmean_2022_2024.mat',model));
        if isfile(f)
            M=load(f,'Meta','testYears','validLinearIndex', ...
                [model '_NRMSE_valid'],[model '_TestSuccess_valid']);
            assert(isequal(double(M.testYears(:))',[2022 2023 2024]) && ...
                isequal(double(M.validLinearIndex(:)),double(P.selectedLinearIndex(:))) && ...
                M.Meta.randomSeed==20260605 && ...
                strcmp(M.Meta.targetDefinition,D.Meta.targetDefinition) && ...
                isequal(size(M.([model '_NRMSE_valid'])),[2000,3]) && ...
                isequal(size(M.([model '_TestSuccess_valid'])),[2000,3]), ...
                'Existing %s output is incompatible; will not reuse it.',model);
            fprintf('%s verified existing output -> skip.\n',model);
            continue;
        end
        fprintf('\n[%d/%d] Training %s, original configuration.\n',j,numel(models),model);
        tm=tic;
        if ismember(model,{'AR','ARIMA','SARIMA'})
            SKLGB_mean_linear(cfg,model);
        else
            feval(['SKLGB_mean_' model],cfg);
        end
        M=load(f,[model '_TestSuccess_valid']);
        ok=M.([model '_TestSuccess_valid']);
        fprintf('MODEL %s: %d/%d successful grid-years; %.1f min.\n', ...
            model,nnz(ok),numel(ok),toc(tm)/60);
        logFile=fullfile(cfg.outputRoot,'training_execution_log.csv');
        isNew=~isfile(logFile);
        fid=fopen(logFile,'a');
        if fid<0,error('Cannot write execution log: %s',logFile);end
        if isNew,fprintf(fid,'Model,SuccessfulGridYears,TotalGridYears,Minutes\n');end
        fprintf(fid,'%s,%d,%d,%.6f\n',model,nnz(ok),numel(ok),toc(tm)/60);
        fclose(fid);
    end
    fprintf('Nine-model training stage complete.\n');
    if strcmp(action,'all')
        SKLGB_mean_analyze(cfg);
    else
        fprintf('Run SKLGB_mean_run(''analyze'') for original-vs-alternative comparisons.\n');
    end
end
end
