function SKLGB_mean_linear(cfg, model)
% Run ORIGINAL AR/ARIMA/SARIMA per-grid estimators on the fixed 2000 cells.
% Model cores are extracted from user-supplied original scripts verbatim:
%   SKLGB_ext_AR_core.m, _ARIMA_core.m, _SARIMA_core.m
% Only the test years, grid subset, I/O and checkpoint orchestration change.
% Numerical model fitting, p selection, test forecasting, failure handling,
% normalization and directional accuracy follow the exact supplied codes.

assert(ismember(model,{'AR','ARIMA','SARIMA'}), 'Unknown linear model.');
assert(isfile(cfg.pilotMasterFile), 'Run preflight to construct sample master.');
D=load(cfg.pilotMasterFile, 'dGWS_valid','originalValidRows', ...
    'validLinearIndex','validLat','validLon','pilotDesignWeight','ym');
Y=double(D.dGWS_valid);
originalValidRows=double(D.originalValidRows(:));
validLinearIndex=double(D.validLinearIndex(:));
validLat=double(D.validLat(:)); validLon=double(D.validLon(:));
pilotDesignWeight=double(D.pilotDesignWeight(:));
ym=double(D.ym(:));
n=size(Y,1); nRound=3; horizon=12;
testYears=[2022 2023 2024];
p_list=[1 6 12 24 36 48 60];
assert(n>=1 && n<=2000 && isequal(size(Y),[n 252]));
assert(numel(ym)==252 && ym(1)==200401 && ym(end)==202412);
assert(numel(unique(validLinearIndex))==n && ...
    numel(originalValidRows)==n && numel(pilotDesignWeight)==n);
assert(all(isfinite(Y(:))), 'Source grid time series contains NaN or Inf.');
assert(cfg.block_num>=1 && cfg.block_num==floor(cfg.block_num));
assert(isequal(double(D.validLinearIndex(:)),validLinearIndex));

sourceNames={'AR','ARIMA','SARIMA'};
sourceHashes={
    '9fbd12d370328ee2edacd18ebc5797912bba3da92a1ed0f0dee55d7de412e0be', ...
    'd451b010f2659a7fb2e887e52e3c5e76b23983f4cdbb1ffe8d81eb80668cdc4d', ...
    'db1ba8729145514bd66fa70c00078b42787088bbd07bcdca98210ed0bf2253d5'};
sourceHash=sourceHashes{strcmp(sourceNames,model)};
sourceFile=fullfile(fileparts(mfilename('fullpath')),...
    ['SKLGB_ext_' model '_core.m']);
assert(isfile(sourceFile),'Missing original extracted model core: %s',sourceFile);
assert(strcmpi(hashFile(sourceFile),sourceHash), ...
    'Source core hash mismatch. Re-extract original model code before running.');

outputDir=fullfile(cfg.outputRoot,model);
if ~isfolder(outputDir),mkdir(outputDir);end
outFile=fullfile(outputDir,sprintf('%s_monthmean_2022_2024.mat',model));
checkpointFile=fullfile(outputDir, ...
    sprintf('%s_monthmean_2022_2024_originalcore_v2_checkpoint.mat',model));
if isfile(outFile)
    error('Final result exists: %s. The launcher must verify or reject it.',outFile);
end
blockSize=ceil(n/cfg.block_num); nBlocks=ceil(n/blockSize);
% SHA-256 of actual data, sampled grid identities, model source and run rules.
runSignature=hashRun(Y,originalValidRows,validLinearIndex, ...
    model,sourceHash,cfg.block_num,testYears,p_list);

if isfile(checkpointFile)
    S=load(checkpointFile);
    assert(isfield(S,'runSignature') && ...
        strcmp(S.runSignature,runSignature), ...
        'Checkpoint conflicts with sample/model code/settings: %s',checkpointFile);
    assert(isequal(size(S.NRMSE),[n,nRound]) && ...
        isequal(size(S.Pred),[n,nRound*horizon]) && ...
        isequal(size(S.done),[nBlocks,1]));
    NRMSE=S.NRMSE; TDA=S.TDA; Pred=S.Pred; ValPred=S.ValPred;
    bestP=S.bestP; ValScore=S.ValScore; Success=S.Success;
    FailureCode=S.FailureCode; FailureID=S.FailureID;
    RMSE=S.RMSE; TestRange=S.TestRange; done=S.done;
    fprintf('Resume source-verified %s: %d/%d blocks.\n',model,nnz(done),nBlocks);
else
    NRMSE=nan(n,3); TDA=nan(n,3); Pred=nan(n,36); ValPred=nan(n,36);
    bestP=nan(n,3); ValScore=nan(n,3); RMSE=nan(n,3);
    TestRange=nan(n,3); Success=false(n,3);
    FailureCode=zeros(n,3,'uint8'); FailureID=cell(n,3);
    done=false(nBlocks,1);
end
if cfg.useParallel
    pool=gcp('nocreate');
    if isempty(pool)
        parpool('local',cfg.numWorkers);
    elseif pool.NumWorkers~=cfg.numWorkers
        error('Existing parallel pool has %d workers, config requires %d.', ...
            pool.NumWorkers,cfg.numWorkers);
    end
end
startTic=tic;
for b=1:nBlocks
    if done(b), fprintf('%s block %d/%d SKIPPED\n',model,b,nBlocks);continue;end
    rows=(b-1)*blockSize+1:min(b*blockSize,n);
    nThis=numel(rows);
    bN=nan(nThis,3);bT=nan(nThis,3);bP=nan(nThis,36);
    bV=nan(nThis,36);bBP=nan(nThis,3);bVS=nan(nThis,3);
    bRMSE=nan(nThis,3);bRange=nan(nThis,3);
    bSuccess=false(nThis,3);bCode=zeros(nThis,3,'uint8');
    bID=cell(nThis,3);
    if cfg.useParallel
        parfor j=1:nThis
            [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
                processOriginalCore(Y(rows(j),:)',model,p_list,testYears,horizon);
            bP(j,:)=pp(:)';bV(j,:)=vp(:)';
            bN(j,:)=nn(:)';bT(j,:)=tt(:)';
            bBP(j,:)=bp(:)';bVS(j,:)=vs(:)';
            bRMSE(j,:)=rr(:)';bRange(j,:)=dr(:)';
            bSuccess(j,:)=ok(:)';bCode(j,:)=fc(:)';bID(j,:)=fid(:)';
        end
    else
        for j=1:nThis
            [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
                processOriginalCore(Y(rows(j),:)',model,p_list,testYears,horizon);
            bP(j,:)=pp(:)';bV(j,:)=vp(:)';
            bN(j,:)=nn(:)';bT(j,:)=tt(:)';
            bBP(j,:)=bp(:)';bVS(j,:)=vs(:)';
            bRMSE(j,:)=rr(:)';bRange(j,:)=dr(:)';
            bSuccess(j,:)=ok(:)';bCode(j,:)=fc(:)';bID(j,:)=fid(:)';
        end
    end
    NRMSE(rows,:)=bN;TDA(rows,:)=bT;Pred(rows,:)=bP;
    ValPred(rows,:)=bV;bestP(rows,:)=bBP;ValScore(rows,:)=bVS;
    RMSE(rows,:)=bRMSE;TestRange(rows,:)=bRange;
    Success(rows,:)=bSuccess;FailureCode(rows,:)=bCode;FailureID(rows,:)=bID;
    done(b)=true;
    temp=[checkpointFile '.tmp.mat'];
    save(temp,'NRMSE','TDA','Pred','ValPred','bestP','ValScore',...
        'RMSE','TestRange','Success','FailureCode','FailureID',...
        'done','runSignature','originalValidRows','-v7.3');
    [okMove,msg]=movefile(temp,checkpointFile,'f');
    if ~okMove,error('Checkpoint replacement failed: %s',msg);end
    fprintf('%s block %d/%d: %d/%d valid grid-years | %.1f min\n', ...
        model,b,nBlocks,nnz(Success(rows,:)),numel(Success(rows,:)),toc(startTic)/60);
end
assert(all(done), 'Not all blocks complete.');

Meta=struct();
Meta.targetDefinition='difference of consecutive 1-degree monthly-mean GWS states';
Meta.model=model;Meta.randomSeed=20260605;Meta.testYears=testYears;
Meta.provenance='ORIGINAL_USER_SUPPLIED_MATLAB_SOURCE';
Meta.sourceCoreSHA256=sourceHash;
Meta.sourceScript=[model '_global_forecast.m'];
Meta.scriptName='SKLGB_mean_linear.m';
Meta.sourceCoreFunction=['SKLGB_ext_' model '_core.m'];
Meta.rolling='Train 2004..T-2; validate T-1; test T; no refit after validation';
Meta.p_list=p_list;Meta.computeTimeMinutes=toc(startTic)/60;
Meta.nValid=n;Meta.pilotOriginalValidRows=originalValidRows;
Meta.nRound=3;Meta.forecastHorizon=12;
Meta.checkpointVersion=2;Meta.runSignature=runSignature;
switch model
    case 'AR'
        Meta.fitMethod='MATLAB ar(Y_train,p,''Ts'',1): default forward-backward';
        Meta.forecastMethod='forecast(bestModel,[train;validation],12)';
    case 'ARIMA'
        Meta.fitMethod='arima(p,0,1); estimate(...,''Display'',''off'')';
        Meta.forecastMethod='forecast(bestModel,12,''Y0'',[train;validation])';
    case 'SARIMA'
        Meta.fitMethod='Original arima(p,0,1)(1,1,1)_12 specification; estimate';
        Meta.forecastMethod='forecast(bestModel,12,''Y0'',[train;validation])';
end

R=struct();R.([model '_NRMSE_valid'])=NRMSE;
R.([model '_TDA_valid'])=TDA;
R.([model '_Pred_valid'])=Pred;
R.([model '_Val_Pred_valid'])=ValPred;
R.([model '_TestSuccess_valid'])=Success;
R.([model '_best_p_valid'])=bestP;
R.([model '_Val_NRMSE_best_valid'])=ValScore;
R.([model '_RMSE_valid'])=RMSE;
R.([model '_TestRange_valid'])=TestRange;
R.([model '_FailureCode_valid'])=FailureCode;
R.([model '_FailureID_valid'])=FailureID;
R.Meta=Meta;R.testYears=testYears(:);
R.validLinearIndex=validLinearIndex;R.validLat=validLat;
R.validLon=validLon;R.originalValidRows=originalValidRows;
R.pilotDesignWeight=pilotDesignWeight;
temp=[outFile '.tmp.mat'];save(temp,'-struct','R','-v7.3');
[okMove,msg]=movefile(temp,outFile,'f');
if ~okMove,error('Final result move failed: %s',msg);end
fprintf('%s ORIGINAL-SOURCE COMPLETED: %d/%d grid-years.\n',...
    model,nnz(Success),numel(Success));
end

function [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
    processOriginalCore(y,model,p_list,testYears,horizon)
switch model
    case 'AR'
        [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
            SKLGB_ext_AR_core(y,p_list,testYears,horizon,1);
    case 'ARIMA'
        [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
            SKLGB_ext_ARIMA_core(y,p_list,0,1,testYears,horizon);
    case 'SARIMA'
        [pp,vp,nn,tt,bp,vs,rr,dr,ok,fc,fid]= ...
            SKLGB_ext_SARIMA_core(y,p_list,0,1,12,1,1,1,testYears,horizon);
    otherwise
        error('Unknown linear model.');
end
end

function h=hashRun(Y,validRows,linearIndex,model,coreHash,blockNum,years,pList)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(int8(uint8([model '_' coreHash])));
md.update(typecast(double(Y(:)),'int8'));
md.update(typecast(double(validRows(:)),'int8'));
md.update(typecast(double(linearIndex(:)),'int8'));
md.update(typecast(double([blockNum years(:)' pList(:)' 2]),'int8'));
b=typecast(int8(md.digest()),'uint8');
h=lower(reshape(dec2hex(b,2).',1,[]));
end

function h=hashFile(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb');
if fid<0,error('Cannot open source core: %s',path);end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    bytes=fread(fid,65536,'*uint8');
    if ~isempty(bytes),md.update(typecast(bytes,'int8'));end
end
b=typecast(int8(md.digest()),'uint8');
h=lower(reshape(dec2hex(b,2).',1,[]));
end
