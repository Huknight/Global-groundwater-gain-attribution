function cfg = SKLGB_mean_config()
% Alternative delta-GWS target: difference of consecutive monthly-mean GWS.
% Only edit cfg.projectRoot and cfg.rawDAFolder below if necessary.
cfg.projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
cfg.rawDAFolder = 'D:\A1_科研数据库\01_原始数据\01_GLDAS\GLDAS-DA CLSM Daily 025D 20030201-to-20251231';

cfg.useParallel = true;
cfg.numWorkers = 15;
cfg.block_num = 20;  % Keep unchanged after training checkpoints exist.
cfg.modelsToRun = {'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};

% Full daily file date coverage is mandatory. Native pixels contribute to
% a monthly state only if their GWS_tavg is valid on EVERY day of that month.
cfg.requireFullDailyCoverage = true;
cfg.checkRawFilesAtPreflight = true;
cfg.checkEveryDailyCoordinates = false;  % First file per month + last day checked.

here = fileparts(mfilename('fullpath'));
cfg.sampleFile = fullfile(here,'SKLGB_hyperparameter_sample2000.mat');
cfg.originalMasterFile = fullfile(cfg.projectRoot, 's1_master_database', ...
    '03_master_database','GLDAS_common_valid_mask_1deg_200401_202412.mat');
cfg.outputRoot = fullfile(cfg.projectRoot,'s3_darft','s1_analysis', ...
    's8_deltaGWS_monthmean_2000grids_singlePC');
cfg.preprocessingDir = fullfile(cfg.outputRoot,'preprocessing');
cfg.nativeMonthDir = fullfile(cfg.preprocessingDir,'native_monthly_mean_states');
cfg.spatialTargetFile = fullfile(cfg.preprocessingDir, ...
    'GLDAS_DA_CLSM10_dGWS_monthlymean_200401_202412.mat');
cfg.pilotMasterFile = fullfile(cfg.outputRoot,'selected2000_monthmean_master.mat');
cfg.nPilot=2000; cfg.nLocal=2000; cfg.machineCount=1; cfg.machineIndex=1;
cfg.rebuildPilotMaster = false;
end
