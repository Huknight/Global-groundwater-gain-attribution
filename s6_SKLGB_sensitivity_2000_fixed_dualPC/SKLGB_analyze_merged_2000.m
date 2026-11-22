function SKLGB_analyze_merged_2000(mergedMat,outputRoot)
% Sample-weighted and original-global-anchored sensitivity comparisons.
% All nine model rankings, six gain classes, and four gain components are
% recomputed for the same 2,000 original cells in every scenario.
if nargin<1 || isempty(mergedMat)
    here=fileparts(mfilename('fullpath'));
    mergedMat=fullfile(here,'merged_results','SKLGB_merged_2000_results.mat');
end
if nargin<2 || isempty(outputRoot),outputRoot=fileparts(mergedMat);end
if ~isfolder(outputRoot),mkdir(outputRoot);end
S=load(mergedMat);
here=fileparts(mfilename('fullpath'));
F=load(fullfile(here,'SKLGB_global_original_reference.mat'));
n=numel(S.originalValidRows);
if n~=2000 || F.populationN~=17284 || ...
        abs(sum(S.designWeight)-F.populationN)>1e-7
    error('Merged sample N / original population reference is inconsistent.');
end
w=double(S.designWeight(:));N=double(F.populationN);
base=double(S.originalAnnualNRMSE9);
if ~isequal(size(base),[2000,9,3]) || any(~isfinite(base(:)))
    error('Original 9-model NRMSE matrix is missing/incomplete.');
end

% Scenarios: unchanged reference, 24 single-model perturbations, and
% 12 coordinated RF/RFX or LSTM/LSTMX perturbations (37 in all).
scNames={'Original_reference'};
scData={base};
modelColumns=[5,8,6,9];
for m=1:4
    for j=2:7
        Q=base;
        Q(:,modelColumns(m),:)=reshape(S.score(:,m,j,:),[n,1,3]);
        scNames{end+1}=[S.models{m} '_only_' S.tagList{m,j}]; %#ok<AGROW>
        scData{end+1}=Q; %#ok<AGROW>
    end
end
for family=1:2
    a=2*family-1; b=2*family;
    for j=2:7
        Q=base;
        Q(:,modelColumns(a),:)=reshape(S.score(:,a,j,:),[n,1,3]);
        Q(:,modelColumns(b),:)=reshape(S.score(:,b,j,:),[n,1,3]);
        scNames{end+1}=['Paired_' S.models{a} '_' S.models{b} '_' S.tagList{a,j}]; %#ok<AGROW>
        scData{end+1}=Q; %#ok<AGROW>
    end
end
if numel(scNames)~=37,error('Expected exactly 37 sensitivity scenarios.');end
allNames={'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};
labels={'2022','2023','2024','Mean_2022_2024'};
classNames={'NoGain','Linear','Nonlinear','Hydroclimatic','MixedL','MixedN'};
componentNames={'Linear','Nonlinear','Hydroclimatic','Mixed'};

% Evaluate the 100% original population reference once against the
% original 2,000-cell sample. It must retain the original six class IDs.
oldQ=cell(4,1);oldClass=cell(4,1);oldWinner=cell(4,1);
oldG=cell(4,1);oldTotal=cell(4,1);
for t=1:4
    if t<4,oldQ{t}=base(:,:,t);else,oldQ{t}=mean(base,3);end
    [oldClass{t},oldWinner{t},oldG{t},oldTotal{t}] = classifyAndDecompose(oldQ{t});
    if t<4,stored=double(S.originalAnnualGainClass(:,t));
    else, stored=double(S.originalMeanGainClass(:));end
    if ~isequal(oldClass{t}(:),stored(:))
        error('Original category reconstruction differs from the frozen sample file.');
    end
end

rankRows=cell(0,9);
classRows=cell(0,7);
shareRows=cell(0,8);
contributionRows=cell(0,8);
for k=1:numel(scNames)
    for t=1:4
        if t<4,Q=scData{k}(:,:,t);else,Q=mean(scData{k},3);end
        if any(~isfinite(Q(:)))
            error('Nonfinite NRMSE in scenario %s, slot %s.',scNames{k},labels{t});
        end
        [C,W,G,Gt]=classifyAndDecompose(Q);
        [~,sampleOrder]=sort(sum(Q.*w,1)/sum(w));
        gModel=zeros(1,9);rankModel=zeros(1,9);
        for m=1:9
            d=sum(w.*(Q(:,m)-oldQ{t}(:,m)))/N;
            gModel(m)=double(F.globalOriginalMeanNRMSE9(m,t))+d;
        end
        [~,globalOrder]=sort(gModel);
        for m=1:9,rankModel(globalOrder(m))=m;end
        weightedMeans=sum(Q.*w,1)/sum(w);
        for m=1:9
            rankRows(end+1,:)={scNames{k},labels{t},allNames{m}, ...
                weightedMeans(m),gModel(m),rankModel(m), ...
                sum(w.*(W==m))/sum(w)*100, ...
                (double(F.globalOriginalWinnerCounts9(m,t)) + ...
                sum(w .* (double(W==m) - double(oldWinner{t}==m)))) / N * 100, ...
                rankModel(m)-find(sampleOrder==m)}; %#ok<AGROW>
        end
        for c=1:6
            sampleCount=sum(w.*double(C==c));
            origCount=sum(w.*double(oldClass{t}==c));
            anchoredCount=double(F.globalOriginalClassCounts6(c,t))+sampleCount-origCount;
            classRows(end+1,:)={scNames{k},labels{t},classNames{c}, ...
                sum(C==c),sum(C==c)/n*100,sampleCount/sum(w)*100, ...
                anchoredCount/N*100}; %#ok<AGROW>
            sumGainC=sum(w(C==c).*Gt(C==c));
            oldClassVector=oldClass{t};
            oldTotalVector=oldTotal{t};
            oldComponentMatrix=oldG{t};
            oldMask=oldClassVector==c;
            sumOldGainC=sum(w(oldMask).*oldTotalVector(oldMask));
            anchoredGainC=double(F.globalOriginalClassGainSums(c,t))+sumGainC-sumOldGainC;
            for d=1:4
                compSum=sum(w(C==c).*G(C==c,d));
                oldComp=sum(w(oldMask).*oldComponentMatrix(oldMask,d));
                anchorComp=double(F.globalOriginalClassComponentSums(c,d,t))+compSum-oldComp;
                if abs(anchoredGainC)>1e-10,share=anchorComp/anchoredGainC*100;
                else,share=NaN;end
                contributionRows(end+1,:)={scNames{k},labels{t}, ...
                    classNames{c},componentNames{d},anchoredGainC, ...
                    anchorComp,share,compSum}; %#ok<AGROW>
            end
        end
        totalAnch=double(F.globalOriginalTotalGainSums(t))+sum(w.*(Gt-oldTotal{t}));
        sComp=zeros(1,4);
        for d=1:4
            sComp(d)=double(F.globalOriginalComponentSums4(d,t)) + ...
                sum(w.*(G(:,d)-oldG{t}(:,d)));
        end
        if abs(sum(sComp)-totalAnch)>1e-6
            error('Gain-component reconstruction failed for %s %s.',scNames{k},labels{t});
        end
        if abs(totalAnch)>1e-10,share=sComp./totalAnch*100;
        else, share=nan(1,4);end
        bestWeighted=sum(w.*min(Q,[],2))/sum(w);
        bestAnch=F.globalOriginalMeanWinnerNRMSE(t) + ...
            sum(w.*(min(Q,[],2)-min(oldQ{t},[],2)))/N;
        shareRows(end+1,:)={scNames{k},labels{t},bestWeighted,bestAnch, ...
            share(1),share(2),share(3),share(4)}; %#ok<AGROW>
    end
end

R=cell2table(rankRows,'VariableNames',{'Scenario','Year','Model', ...
    'SampleWeightedMeanNRMSE','GlobalAnchoredMeanNRMSE', ...
    'GlobalAnchoredRank','SampleWeightedWinnerPct', ...
    'GlobalAnchoredWinnerPct','RankVsWeightedSample'});
% Guard against integer division silently rounding winner percentages to zero.
for iScenario = 1:numel(scNames)
    for iYear = 1:numel(labels)
        ix = (iScenario - 1) * numel(labels) * 9 + (iYear - 1) * 9 + (1:9);
        if abs(sum(R.GlobalAnchoredWinnerPct(ix)) - 100) > 1e-8
            error('Anchored winner percentages do not sum to 100 for %s %s.', ...
                scNames{iScenario}, labels{iYear});
        end
    end
end
writetable(R,fullfile(outputRoot,'SKLGB_model_rankings_2000.csv'));
R=cell2table(classRows,'VariableNames',{'Scenario','Year','GainClass', ...
    'RawSampleCount','RawSamplePct','SampleWeightedPct', ...
    'GlobalAnchoredPct'});
writetable(R,fullfile(outputRoot,'SKLGB_gain_class_comparison_2000.csv'));
R=cell2table(shareRows,'VariableNames',{'Scenario','Year', ...
    'SampleWeightedBestNRMSE','GlobalAnchoredBestNRMSE', ...
    'LinearSharePct','NonlinearSharePct','HydroSharePct','MixedSharePct'});
writetable(R,fullfile(outputRoot,'SKLGB_global_component_shares_2000.csv'));
R=cell2table(contributionRows,'VariableNames',{'Scenario','Year', ...
    'GainClass','Component','GlobalAnchoredClassTotalGainSum', ...
    'GlobalAnchoredComponentSum','GlobalAnchoredWithinClassPct', ...
    'SampleWeightedComponentSum'});
writetable(R,fullfile(outputRoot,'SKLGB_gain_components_by_class_2000.csv'));

fprintf('ANALYSIS COMPLETE: 37 scenarios x 4 time slices.\n');
fprintf('Model ranks, six gain classes, and gain components are saved in %s\n',outputRoot);
end

function [C,W,G,total]=classifyAndDecompose(Q)
% Faithful to the six-category algorithm in SKLGB_pilot10_summary.m and
% G-component telescoping relations in s2_gain_attribution3.m.
n=size(Q,1);
if size(Q,2)~=9 || any(~isfinite(Q(:)))
    error('classifyAndDecompose requires a finite N-by-9 NRMSE matrix.');
end
[~,W]=min(Q,[],2);
[L,li]=min(Q(:,1:3),[],2);
[N,ni]=min(Q(:,4:6),[],2);
[X,xi]=min(Q(:,7:9),[],2);
[~,groups]=sort([L N X],2,'ascend');
C=zeros(n,1);
best=groups(:,1);second=groups(:,2);
C(best==1 & li==1)=1;
C(best==1 & li~=1)=2;
C(best==2)=3;
C(best==3 & second==2 & xi==ni)=4;
C(best==3 & second==1)=5;
C(best==3 & second==2 & xi~=ni)=6;
if any(C==0),error('Unclassified cells found in gain decomposition.');end
AR=Q(:,1);B=min(Q,[],2);
total=AR-B;
G=zeros(n,4);
t=C==2;G(t,1)=AR(t)-B(t);
t=C==3;G(t,1)=AR(t)-L(t);G(t,2)=L(t)-N(t);
t=C==4;G(t,1)=AR(t)-L(t);G(t,2)=L(t)-N(t);G(t,3)=N(t)-X(t);
t=C==5;G(t,1)=AR(t)-L(t);G(t,4)=L(t)-X(t);
t=C==6;G(t,1)=AR(t)-L(t);G(t,2)=L(t)-N(t);G(t,4)=N(t)-X(t);
if max(abs(sum(G,2)-total))>1e-8
    error('Non-telescoping gain components found in sensitivity analysis.');
end
end
