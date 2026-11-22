function SKLGB_pilot10_run(machineIndex, action)
% FIXED 2,000-grid / two-computer sensitivity experiment.
% machineIndex MUST be 1 (PC1) or 2 (PC2); no sample-size switch.
% action='preflight' verifies allocation/plots without training;
% action='train' runs 24 OAT configurations per computer.
if nargin < 1 || ~isscalar(machineIndex) || ~ismember(machineIndex,[1,2])
    error('Use SKLGB_PC1_train or SKLGB_PC2_train (PC number must be 1 or 2).');
end
if nargin < 2
    action = 'train';
end
if ~ismember(action, {'preflight','train'})
    error('Action must be preflight or train.');
end
% User-editable settings are kept in this one section only.
% Required: MATLAB Statistics and Machine Learning Toolbox, Deep Learning
% Toolbox, Parallel Computing Toolbox if useParallel=true.

%% ======== EDIT ONLY LOCAL COMPUTER SETTINGS ========
projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';  % edit if path differs
sampleFile = fullfile(fileparts(mfilename('fullpath')), ...
    'SKLGB_hyperparameter_sample2000.mat');
nPilot = 2000;                % LOCKED: always 2,000 selected grid cells
machineCount = 2;             % LOCKED: PC1 and PC2, 1,000 each
useParallel = true;
numWorkers = 15;              % adjust for computer resources
block_num = 30;               % LOCKED once a checkpoint exists
mode = 'perturb_only';        % LOCKED: no original baseline refitting
modelsToRun = {'RF','RFX','LSTM','LSTMX'};
rebuildPilotMaster = false;
plotSampleMap = true;
%% =========================================

cfg = struct();
cfg.projectRoot = projectRoot;
cfg.originalMasterFile = fullfile(projectRoot, 's1_master_database', ...
    '03_master_database', 'GLDAS_common_valid_mask_1deg_200401_202412.mat');
cfg.sampleFile = sampleFile;
cfg.outputRoot = fullfile(projectRoot, 's3_darft', 's1_analysis', ...
    sprintf('s6_hyperparameter_sensitivity_%dgrids_machine%dof%d', ...
    nPilot,machineIndex,machineCount));
cfg.pilotMasterFile = fullfile(cfg.outputRoot,'pilot10_master.mat');
if machineCount < 1 || machineCount ~= floor(machineCount) || ...
        machineIndex < 1 || machineIndex > machineCount || ...
        machineIndex ~= floor(machineIndex) || machineCount > nPilot
    error('Invalid machineCount or machineIndex.');
end
cfg.nPilot = nPilot;
cfg.machineCount = machineCount;
cfg.machineIndex = machineIndex;
cfg.nLocal = floor((nPilot-machineIndex)/machineCount)+1;
cfg.plotSampleMap = plotSampleMap;
cfg.useParallel = useParallel;
cfg.numWorkers = numWorkers;
cfg.block_num = block_num;
cfg.rebuildPilotMaster = rebuildPilotMaster;
if ~ismember(mode,{'all','baseline_only','perturb_only'})
    error('mode must be all, baseline_only or perturb_only.');
end
validModels = {'RF','RFX','LSTM','LSTMX'};
if ~all(ismember(modelsToRun,validModels)) || isempty(modelsToRun)
    error('modelsToRun must be a nonempty subset of RF,RFX,LSTM,LSTMX.');
end
if ~isfolder(cfg.outputRoot), mkdir(cfg.outputRoot); end
if nPilot~=2000 || machineCount~=2 || ~strcmp(mode,'perturb_only')
    error('This release is locked to 2,000 grids, two PCs, perturb_only.');
end
fprintf('\n==============================================================\n');
fprintf('SKLGB FIXED 2000-grid run: %d globally, PC %d/%d with %d cells\n', ...
    cfg.nPilot,cfg.machineIndex,cfg.machineCount,cfg.nLocal);
fprintf('Mode: %s | Parallel workers requested: %d\n', mode, cfg.numWorkers);
fprintf('All parameter variants reselect p (and RFX x) by validation NRMSE.\n');
fprintf('Experiment seed: ORIGINAL run 20260605; original point IDs retained.\n');
fprintf('==============================================================\n');
SKLGB_validate_selection(cfg.sampleFile);
SKLGB_pilot10_prepare_master(cfg);
if cfg.plotSampleMap
    SKLGB_plot_sampling_map(cfg);
end
if strcmp(action,'preflight')
    fprintf('PREFLIGHT PASSED: PC %d/2, %d assigned original grid cells. No training done.\n', ...
        cfg.machineIndex,cfg.nLocal);
    return;
end

% Every model runs ORIGINAL baseline first; all six OAT deviations follow.
% The seven variants per model share the original p/x search ranges.
runConfigs = struct('model',{},'tag',{},'trees',{},'leaf',{},'hidden',{},'lr',{});
for m = 1:numel(modelsToRun)
    model = modelsToRun{m};
    if ismember(model,{'RF','RFX'})
        vals = { 'baseline',50,5; 'trees_025',25,5; ...
            'trees_100',100,5; 'trees_200',200,5; ...
            'leaf_001',50,1; 'leaf_010',50,10; 'leaf_020',50,20 };
        for j=1:size(vals,1)
            runConfigs(end+1) = struct('model',model,'tag',vals{j,1}, ...
                'trees',vals{j,2},'leaf',vals{j,3},'hidden',32,'lr',0.01); %#ok<AGROW>
        end
    else
        vals = { 'baseline',32,0.01; 'hidden_016',16,0.01; ...
            'hidden_064',64,0.01; 'hidden_128',128,0.01; ...
            'lr_0p001',32,0.001; 'lr_0p003',32,0.003; ...
            'lr_0p030',32,0.03 };
        for j=1:size(vals,1)
            runConfigs(end+1) = struct('model',model,'tag',vals{j,1}, ...
                'trees',50,'leaf',5,'hidden',vals{j,2},'lr',vals{j,3}); %#ok<AGROW>
        end
    end
end

% First complete all model baselines before any perturbation.
runOrder = [find(strcmp({runConfigs.tag}, 'baseline')), ...
    find(~strcmp({runConfigs.tag}, 'baseline'))];
if strcmp(mode,'baseline_only')
    runOrder = find(strcmp({runConfigs.tag}, 'baseline'));
elseif strcmp(mode,'perturb_only')
    runOrder = find(~strcmp({runConfigs.tag}, 'baseline'));
end

for seq = 1:numel(runOrder)
    k = runOrder(seq);
    runCfg = cfg;
    runCfg.model = runConfigs(k).model;
    runCfg.tag = runConfigs(k).tag;
    runCfg.trees = runConfigs(k).trees;
    runCfg.leaf = runConfigs(k).leaf;
    runCfg.hidden = runConfigs(k).hidden;
    runCfg.lr = runConfigs(k).lr;
    outDir = fullfile(cfg.outputRoot,runCfg.model,runCfg.tag);
    if ~isfolder(outDir), mkdir(outDir); end
    modelOutputFile = SKLGB_pilot10_result_path(cfg.outputRoot,runCfg.model,runCfg.tag);
    fprintf('\n[%d/%d] %s | %s | trees=%d leaf=%d hidden=%d LR=%g\n', ...
        seq,numel(runOrder),runCfg.model,runCfg.tag,runCfg.trees, ...
        runCfg.leaf,runCfg.hidden,runCfg.lr);
    taskTic = tic;
    taskStatus = 'executed';
    if isfile(modelOutputFile)
        taskStatus = 'skipped_existing';
        A = load(modelOutputFile,'Meta','validLinearIndex');
        M = load(cfg.pilotMasterFile,'originalLinearIndex','originalValidRows');
        if ~isfield(A,'Meta') || ~isfield(A.Meta,'pilotTag') || ...
                ~strcmp(A.Meta.pilotTag,runCfg.tag) || ...
                ~isequal(double(A.validLinearIndex(:)),double(M.originalLinearIndex(:))) || ...
                ~isequal(double(A.Meta.pilotOriginalValidRows(:)),double(M.originalValidRows(:)))
            error('Existing result is incompatible. Inspect: %s', modelOutputFile);
        end
        if ~strcmp(A.Meta.experimentName,'original') || ...
                A.Meta.randomSeed ~= 20260605 || A.Meta.nValid ~= cfg.nLocal
            error('Existing run identity mismatch: %s',modelOutputFile);
        end
        if ismember(runCfg.model,{'RF','RFX'})
            if A.Meta.rfParam.NumLearningCycles ~= runCfg.trees || ...
                    A.Meta.rfParam.MinLeafSize ~= runCfg.leaf
                error('Existing RF hyperparameters mismatch: %s',modelOutputFile);
            end
        else
            if A.Meta.lstmParam.numHiddenUnits ~= runCfg.hidden || ...
                    A.Meta.lstmParam.initialLearnRate ~= runCfg.lr
                error('Existing LSTM hyperparameters mismatch: %s',modelOutputFile);
            end
        end
        fprintf('Existing compatible final output: skipping training.\n');
    else
        tStart = tic;
        switch runCfg.model
            case 'RF', SKLGB_pilot_RF(runCfg);
            case 'RFX', SKLGB_pilot_RFX(runCfg);
            case 'LSTM', SKLGB_pilot_LSTM(runCfg);
            case 'LSTMX', SKLGB_pilot_LSTMX(runCfg);
        end
        fprintf('Configuration wall-time (including final save): %.3f min\n', ...
            toc(tStart)/60);
    end
    % Persist wall time and model fit time after every completed configuration.
    Mdone = load(modelOutputFile,'Meta', [runCfg.model '_TestSuccess_valid']);
    succField = [runCfg.model '_TestSuccess_valid'];
    nSuccess = nnz(Mdone.(succField));
    wallLog = fullfile(cfg.outputRoot,'pilot10_execution_log.csv');
    isNewLog = ~isfile(wallLog);
    fid = fopen(wallLog,'a');
    if fid<0, error('Cannot append wall-time log: %s',wallLog); end
    if isNewLog
        fprintf(fid,'Model,Configuration,Status,WallTimeMinutes,MetaComputeTimeMinutes,SuccessfulGridYears,TotalGridYears\n');
    end
    fprintf(fid,'%s,%s,%s,%.6f,%.6f,%d,%d\n', ...
        runCfg.model,runCfg.tag,taskStatus,toc(taskTic)/60, ...
        Mdone.Meta.computeTimeMinutes,nSuccess,cfg.nLocal*3);
    fclose(fid);
    % Write partial summaries even when some later configurations have not run.
    SKLGB_pilot10_summary(cfg);
end

SKLGB_export_compact_results(cfg.outputRoot);
Sdone=load(fullfile(cfg.outputRoot,'compact_sensitivity_results.mat'), ...
    'wasFit','successCount');
if ~all(all(Sdone.wasFit(:,2:7))) || ...
        any(any(Sdone.successCount(:,2:7)~=cfg.nLocal*3))
    error('Not all 24 perturbation outputs are complete. Check the log.');
end
fprintf('ALL 24 PERTURBATION CONFIGURATIONS COMPLETE for PC %d.\n', ...
    cfg.machineIndex);
fprintf('\nSample run completed. Summary folder: %s\n',cfg.outputRoot);
fprintf('Main files: pilot10_model_metrics.csv, pilot10_attribution.csv, ...\n');
fprintf('            pilot10_baseline_reproduction.csv\n');
end

function path = SKLGB_pilot10_result_path(root,model,tag)
switch model
    case 'RF', name='RF_univariate_global_200401_202412.mat';
    case 'RFX', name='RFX_multivariate_global_200401_202412.mat';
    case 'LSTM', name='LSTM_univariate_global_200401_202412.mat';
    case 'LSTMX', name='LSTMX_multivariate_global_200401_202412.mat';
    otherwise, error('Unknown model: %s',model);
end
path=fullfile(root,model,tag,name);
end
