function cfg = SKLGB_ext_config()
% One editable settings file for the 2000-grid / one-workstation experiment.
%% ===== CHANGE PROJECT ROOT IF DIFFERENT =====
cfg.projectRoot = 'D:\P1_SKLGB\M_manuscript1_revised';
cfg.useParallel = true;
cfg.numWorkers = 15;
cfg.block_num = 20;  % DO NOT change after checkpoint is written
% All nine models preserve the actual user-supplied original MATLAB cores.
% One PC, one parallel pool, models executed sequentially (no two-PC split).
% Place linear models first; existing completed nonlinear models are skipped.
cfg.modelsToRun = {'AR','ARIMA','SARIMA','RF','RFX','LSTM','LSTMX','SVM','SVMX'};
%% ===== NO OTHER EDITS NEEDED =====
here = fileparts(mfilename('fullpath'));
cfg.sampleFile = fullfile(here,'SKLGB_hyperparameter_sample2000.mat');
cfg.originalMasterFile = fullfile(cfg.projectRoot,'s1_master_database', ...
    '03_master_database','GLDAS_common_valid_mask_1deg_200401_202412.mat');
cfg.outputRoot = fullfile(cfg.projectRoot,'s3_darft','s1_analysis', ...
    's7_extended_2019_2021_2000grids_singlePC');
cfg.pilotMasterFile = fullfile(cfg.outputRoot,'selected2000_master.mat');
cfg.nPilot = 2000;
cfg.nLocal = 2000;
cfg.machineCount = 1;
cfg.machineIndex = 1;
cfg.rebuildPilotMaster = false;
end
