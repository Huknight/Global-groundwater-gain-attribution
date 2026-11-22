function SKLGB_pilot10_summary(cfg)
% Summarize all currently finished baseline/OAT configurations.
% Pilot data are diagnostic only: ten grid cells do NOT estimate a global CI.
Q = load(cfg.pilotMasterFile,'pilotBaselineNRMSE9','pilotMeanGainClass', ...
    'pilotAnnualGainClass','originalValidRows');
base = double(Q.pilotBaselineNRMSE9);
n = size(base,1);
if ~isequal(size(base),[n,9,3])
    error('Expected baseline NRMSE array [pilot cell, nine models, three years].');
end
modelNames={'RF','RFX','LSTM','LSTMX'};
modelIndices=[5,8,6,9];
tagsRF={'baseline','trees_025','trees_100','trees_200', ...
    'leaf_001','leaf_010','leaf_020'};
tagsLSTM={'baseline','hidden_016','hidden_064','hidden_128', ...
    'lr_0p001','lr_0p003','lr_0p030'};

modelMetricCells=cell(0,11);
baselineCells=cell(0,7);
modelResults=struct();
for m=1:numel(modelNames)
    model=modelNames{m};
    col=modelIndices(m);
    if m<=2, tags=tagsRF; else, tags=tagsLSTM; end
    for k=1:numel(tags)
        tag=tags{k};
        path=resultFile(cfg.outputRoot,model,tag);
        if ~isfile(path), continue; end
        V = load(path,[model '_NRMSE_valid'],[model '_best_p_valid'],'Meta');
        mat = double(V.([model '_NRMSE_valid']));
        P = double(V.([model '_best_p_valid']));
        if ~isequal(size(mat),[n,3]) || ~isequal(size(P),[n,3])
            error('Invalid pilot model output dimensions: %s',path);
        end
        if ~isfield(V,'Meta') || ~strcmp(V.Meta.experimentName,'original') || ...
                ~isfield(V.Meta,'pilotTag') || ~strcmp(V.Meta.pilotTag,tag)
            error('Incorrect experiment identity in: %s',path);
        end
        if ~isequal(double(V.Meta.pilotOriginalValidRows(:)),double(Q.originalValidRows(:)))
            error('Grid identity mismatch in: %s',path);
        end
        key=matlab.lang.makeValidName([model '_' tag]);
        modelResults.(key)=mat;
        if strcmp(model,'RFX')
            X = load(path,'RFX_best_x_valid');
            xChosen=double(X.RFX_best_x_valid);
        else
            xChosen=nan(size(P));
        end
        ref = squeeze(base(:,col,:));
        if strcmp(tag,'baseline')
            dif = mat - ref;
            ok = isfinite(dif);
            if any(ok(:))
                maxAbsDiff=max(abs(dif(ok)));
                meanAbsDiff=mean(abs(dif(ok)));
            else
                maxAbsDiff=NaN; meanAbsDiff=NaN;
            end
            baselineCells(end+1,:)={model,nnz(ok),n*3,maxAbsDiff,meanAbsDiff, ...
                nnz(~isfinite(mat)), V.Meta.randomSeed}; %#ok<AGROW>
        end
        for yr=0:3
            if yr==0
                v = mat(:);
                r = ref(:);
                pp=P(:);
                xx=xChosen(:);
                yearLabel='2022-2024';
            else
                v=mat(:,yr);
                r=ref(:,yr);
                pp=P(:,yr);
                xx=xChosen(:,yr);
                yearLabel=num2str(2021+yr);
            end
            finite=isfinite(v) & isfinite(r);
            if any(finite)
                meanVal=mean(v(finite));
                meanRef=mean(r(finite));
                meanDelta=mean(v(finite)-r(finite));
            else
                meanVal=NaN;meanRef=NaN;meanDelta=NaN;
            end
            meanP=mean(pp(isfinite(pp)));
            if strcmp(model,'RFX')
                meanX=mean(xx(isfinite(xx)));
            else
                meanX=NaN;
            end
            modelMetricCells(end+1,:)={model,tag,yearLabel,nnz(finite), ...
                numel(v),meanRef,meanVal,meanDelta,meanP,meanX, ...
                V.Meta.computeTimeMinutes}; %#ok<AGROW>
        end
    end
end
metricNames={'Model','Configuration','Year','SuccessCount','TotalCount', ...
    'OriginalMeanNRMSE','PilotMeanNRMSE','PairedMeanDeltaNRMSE', ...
    'MeanSelectedP','MeanSelectedX','ComputeTimeMinutes'};
if ~isempty(modelMetricCells)
    T=cell2table(modelMetricCells,'VariableNames',metricNames);
    writetable(T,fullfile(cfg.outputRoot,'pilot10_model_metrics.csv'));
    fprintf('\nPILOT METRICS: %d model-configuration-year rows saved.\n',height(T));
    disp(T(strcmp(T.Year,'2022-2024'),:));
end
if ~isempty(baselineCells)
    T=cell2table(baselineCells,'VariableNames', ...
        {'Model','ComparedValues','TotalValues','MaxAbsNRMSEDifference', ...
        'MeanAbsNRMSEDifference','FailedGridYears','RandomSeed'});
    writetable(T,fullfile(cfg.outputRoot,'pilot10_baseline_reproduction.csv'));
    fprintf('\nORIGINAL BASELINE REPRODUCTION CHECK:\n');
    disp(T);
end

% Recompute winner-conditioned attribution for (a) original ten-cell
% baseline, (b) single-model deviations, and (c) matched paired deviations.
scenarioLabels={'Original_reference'};
scenarioScores={base};
for m=1:numel(modelNames)
    model=modelNames{m}; col=modelIndices(m);
    if m<=2,tags=tagsRF;else,tags=tagsLSTM;end
    for k=2:numel(tags)
        tag=tags{k}; key=matlab.lang.makeValidName([model '_' tag]);
        if isfield(modelResults,key)
            newScore=base;
            newScore(:,col,:)=reshape(modelResults.(key),[n,1,3]);
            scenarioLabels{end+1}=[model '_only_' tag]; %#ok<AGROW>
            scenarioScores{end+1}=newScore; %#ok<AGROW>
        end
    end
end
for k=2:numel(tagsRF)
    tag=tagsRF{k};
    key1=matlab.lang.makeValidName(['RF_' tag]);
    key2=matlab.lang.makeValidName(['RFX_' tag]);
    if isfield(modelResults,key1) && isfield(modelResults,key2)
        newScore=base;
        newScore(:,5,:)=reshape(modelResults.(key1),[n,1,3]);
        newScore(:,8,:)=reshape(modelResults.(key2),[n,1,3]);
        scenarioLabels{end+1}=['Paired_RF_RFX_' tag]; %#ok<AGROW>
        scenarioScores{end+1}=newScore; %#ok<AGROW>
    end
end
for k=2:numel(tagsLSTM)
    tag=tagsLSTM{k};
    key1=matlab.lang.makeValidName(['LSTM_' tag]);
    key2=matlab.lang.makeValidName(['LSTMX_' tag]);
    if isfield(modelResults,key1) && isfield(modelResults,key2)
        newScore=base;
        newScore(:,6,:)=reshape(modelResults.(key1),[n,1,3]);
        newScore(:,9,:)=reshape(modelResults.(key2),[n,1,3]);
        scenarioLabels{end+1}=['Paired_LSTM_LSTMX_' tag]; %#ok<AGROW>
        scenarioScores{end+1}=newScore; %#ok<AGROW>
    end
end

atrRows=cell(0,17);
for k=1:numel(scenarioScores)
    dat=scenarioScores{k};
    for yr=0:3
        if yr==0
            Q9=mean(dat,3); % all three years complete for required 9 models
            yearLabel='Mean_2022_2024';
        else
            Q9=dat(:,:,yr);
            yearLabel=num2str(2021+yr);
        end
        mask=all(isfinite(Q9),2);
        Q9=Q9(mask,:);
        nn=size(Q9,1);
        if nn>0
            [class,winner]=classifyGain(Q9);
            % Column differences: RF-RFX and LSTM-LSTMX; positive means
            % the augmented (X) counterpart has a lower prediction error.
            categories=arrayfun(@(x)nnz(class==x),1:6);
            winnerMeans=mean(min(Q9,[],2));
            rfMean=mean(Q9(:,5)); rfxMean=mean(Q9(:,8));
            lstmMean=mean(Q9(:,6)); lstmxMean=mean(Q9(:,9));
            rfGain=mean(Q9(:,5)-Q9(:,8));
            lstmGain=mean(Q9(:,6)-Q9(:,9));
            if yr==0
                originalQ=mean(base,3);
            else
                originalQ=base(:,:,yr);
            end
            [oldClass,oldWinner]=classifyGain(originalQ(mask,:));
            winnerRetention=mean(winner==oldWinner);
            classRetention=mean(class==oldClass);
        else
            categories=zeros(1,6);
            winnerMeans=NaN;
            rfMean=NaN;rfxMean=NaN;lstmMean=NaN;lstmxMean=NaN;
            rfGain=NaN;lstmGain=NaN;
            winnerRetention=NaN;classRetention=NaN;
        end
        atrRows(end+1,:)={scenarioLabels{k},yearLabel,nn,n, ...
            categories(1),categories(2),categories(3),categories(4), ...
            categories(5),categories(6),winnerMeans,rfMean,rfxMean, ...
            lstmMean,lstmxMean,rfGain,lstmGain}; %#ok<AGROW>
        retentionRows{size(atrRows,1),1}=scenarioLabels{k}; %#ok<AGROW>
        retentionRows{size(atrRows,1),2}=yearLabel;
        retentionRows{size(atrRows,1),3}=winnerRetention;
        retentionRows{size(atrRows,1),4}=classRetention;
    end
end
if ~isempty(atrRows)
    A=cell2table(atrRows,'VariableNames',{'Scenario','Year','CompleteCells', ...
        'TotalCells','NoGain','Linear','Nonlinear','Hydroclimatic', ...
        'MixedL','MixedN','MeanWinnerNRMSE','MeanRF','MeanRFX', ...
        'MeanLSTM','MeanLSTMX','MeanRFminusRFX','MeanLSTMminusLSTMX'});
    writetable(A,fullfile(cfg.outputRoot,'pilot10_attribution.csv'));
    R=cell2table(retentionRows,'VariableNames',{'Scenario','Year', ...
        'WinnerRetention','GainClassRetention'});
    writetable(R,fullfile(cfg.outputRoot,'pilot10_attribution_retention.csv'));
    fprintf('PILOT ATTRIBUTION: %d scenario-year rows saved.\n',height(A));
    if numel(scenarioScores)>1
        show=strcmp(A.Year,'Mean_2022_2024');
        disp(A(show,{'Scenario','CompleteCells','NoGain','Linear', ...
            'Nonlinear','Hydroclimatic','MixedL','MixedN','MeanWinnerNRMSE'}));
    end
end

% Assert the classification implementation reproduces the original stored
% annual and 3-year mean class codes, not just the aggregated class counts.
[referenceMeanClass,~]=classifyGain(mean(base,3));
if ~isequal(referenceMeanClass(:),double(Q.pilotMeanGainClass(:)))
    error('Mean gain-attribution reconstruction differs from the sample MAT.');
end
for yr=1:3
    [referenceYearClass,~]=classifyGain(base(:,:,yr));
    if ~isequal(referenceYearClass(:),double(Q.pilotAnnualGainClass(:,yr)))
        error('Annual %d gain-attribution reconstruction differs from sample MAT.',2021+yr);
    end
end
end

function path=resultFile(root,model,tag)
switch model
    case 'RF', name='RF_univariate_global_200401_202412.mat';
    case 'RFX', name='RFX_multivariate_global_200401_202412.mat';
    case 'LSTM', name='LSTM_univariate_global_200401_202412.mat';
    case 'LSTMX', name='LSTMX_multivariate_global_200401_202412.mat';
    otherwise, error('Unknown model');
end
path=fullfile(root,model,tag,name);
end

function [classes,winner]=classifyGain(scores)
% Exact original six-class winner-conditioned attribution conditions.
% 1 No gain / 2 Linear / 3 Nonlinear / 4 Hydroclimatic / 5 Mixed-L / 6 Mixed-N
n=size(scores,1);
classes=nan(n,1); winner=nan(n,1);
if n==0,return;end
if size(scores,2)~=9 || any(~isfinite(scores(:)))
    error('classifyGain requires finite n-by-9 NRMSE values.');
end
[~,winner]=min(scores,[],2);
[linear,linearIdx]=min(scores(:,1:3),[],2);
[nonlinear,nonlinearIdx]=min(scores(:,4:6),[],2);
[augmented,augIdx]=min(scores(:,7:9),[],2);
groupScore=[linear,nonlinear,augmented];
[~,groupOrder]=sort(groupScore,2,'ascend');
for i=1:n
    best=groupOrder(i,1);
    second=groupOrder(i,2);
    if best==1
        if linearIdx(i)==1,classes(i)=1;else,classes(i)=2;end
    elseif best==2
        classes(i)=3;
    elseif best==3
        if second==2 && augIdx(i)==nonlinearIdx(i)
            classes(i)=4;
        elseif second==1
            classes(i)=5;
        else
            classes(i)=6;
        end
    end
end
end

