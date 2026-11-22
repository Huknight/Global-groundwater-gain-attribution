function SKLGB_ext_analyze(cfg)
% Compare six separate annual test years on the SAME fixed 2000 cells.
% No six-year winner, no population-anchored extrapolation.
if nargin<1, cfg=SKLGB_ext_config(); end
P=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex', ...
    'selectedLat','selectedLon','designWeight','annualNRMSE9', ...
    'annualGainClass','annualWinnerModel');
models={'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};
years=2019:2024; n=2000;
assert(isequal(size(P.annualNRMSE9),[n,9,3]));
weights=double(P.designWeight(:));
assert(all(weights>0) && abs(sum(weights)-17284)<1e-6);
E=nan(n,9,6);
E(:,:,4:6)=double(P.annualNRMSE9);
provenance=cell(9,1);
for m=1:9
    name=models{m};
    path=fullfile(cfg.outputRoot,name, ...
        sprintf('%s_extended_2019_2021.mat',name));
    if ~isfile(path), error('Missing result for %s: %s',name,path); end
    A=load(path,[name '_NRMSE_valid'],[name '_TestSuccess_valid'], ...
        'validLinearIndex','testYears','Meta');
    Z=double(A.([name '_NRMSE_valid']));
    OK=A.([name '_TestSuccess_valid']);
    assert(isequal(size(Z),[n,3]) && isequal(size(OK),[n,3]));
    assert(isequal(double(A.validLinearIndex(:)),double(P.selectedLinearIndex(:))));
    assert(isequal(double(A.testYears(:))',[2019 2020 2021]));
    assert(A.Meta.randomSeed==20260605,'Seed mismatch for %s',name);
    Z(~OK)=NaN;
    E(:,m,1:3)=reshape(Z,[n,1,3]);
    provenance{m}=A.Meta.scriptName;
end
if any(~isfinite(E(:)))
    common = squeeze(all(isfinite(E),2));
    fprintf('WARNING: nonfinite model outputs. Common nine-model coverage per year: ');
    fprintf('%d ',sum(common,1)); fprintf('\n');
    if any(sum(common,1)<0.95*n)
        error('At least one year has <95%% complete common cells; inspect failures.');
    end
else
    common=true(n,6);
end
% Verify exact attribution implementation by reproducing existing 2022-24.
for j=1:3
    [win,cat,~,~]=attributeYear(E(:,:,j+3));
    good=common(:,j+3);
    assert(all(win(good)==double(P.annualWinnerModel(good,j))));
    assert(all(cat(good)==double(P.annualGainClass(good,j))));
end
fprintf('2022-2024 original sample: reproduced exact winner and gain category.\n');

modelRows=cell(0,6); winnerRows=cell(0,5);
classRows=cell(0,5); componentRows=cell(0,7);
overallRows=cell(0,5); summaryRows=cell(0,13);
Winner=nan(n,6); Class=nan(n,6); Components=nan(n,4,6);
Period = cell(6,1);
for iy=1:6
    yr=years(iy); X=E(:,:,iy); mask=common(:,iy);
    Period{iy}='2019-2021';
    if yr>=2022, Period{iy}='2022-2024'; end
    [win,cls,comp,gain]=attributeYear(X);
    win(~mask)=NaN; cls(~mask)=NaN; comp(~mask,:)=NaN;
    Winner(:,iy)=win; Class(:,iy)=cls; Components(:,:,iy)=comp;
    wi=weights(mask);
    scores=nan(1,9);
    for m=1:9, scores(m)=sum(weights(mask).*X(mask,m))/sum(wi); end
    [~,ii]=sort(scores,'ascend');
    for m=1:9
        rank=find(ii==m,1);
        winnerPct=100*sum(weights(mask & win==m))/sum(wi);
        modelRows(end+1,:)={yr,Period{iy},models{m},scores(m),rank,sum(mask)}; %#ok<AGROW>
        winnerRows(end+1,:)={yr,Period{iy},models{m},winnerPct, ...
            sum(mask & win==m)}; %#ok<AGROW>
    end
    for c=1:6
        cc=mask & cls==c;
        p=100*sum(weights(cc))/sum(wi);
        classRows(end+1,:)={yr,Period{iy},categoryName(c),p,nnz(cc)}; %#ok<AGROW>
        if c==1,continue;end
        total=sum(weights(cc).*gain(cc));
        for k=1:4
            if total>0, share=100*sum(weights(cc).*comp(cc,k))/total;
            else, share=NaN; end
            componentRows(end+1,:)={yr,Period{iy},categoryName(c), ...
                componentName(k),share,sum(weights(cc).*comp(cc,k)), ...
                total}; %#ok<AGROW>
        end
    end
    tot=sum(weights(mask).*gain(mask));
    overall=nan(1,4);
    for k=1:4
        if tot>0, overall(k)=100*sum(weights(mask).*comp(mask,k))/tot; end
        overallRows(end+1,:)={yr,Period{iy},componentName(k), ...
            overall(k),sum(weights(mask).*comp(mask,k))}; %#ok<AGROW>
    end
    hydro=mask & cls==4;
    hydroTot=sum(weights(hydro).*gain(hydro));
    if hydroTot>0
        hydroHistory=100*sum(weights(hydro).*(comp(hydro,1)+comp(hydro,2)))/hydroTot;
        hydroFinal=100*sum(weights(hydro).*comp(hydro,3))/hydroTot;
    else
        hydroHistory=NaN; hydroFinal=NaN;
    end
    structPct=100*sum(weights(mask & (cls==2 | cls==3)))/sum(wi);
    hydroPct=100*sum(weights(hydro))/sum(wi);
    mixedPct=100*sum(weights(mask & (cls==5 | cls==6)))/sum(wi);
    bestNRMSE=sum(weights(mask).*min(X(mask,:),[],2))/sum(wi);
    summaryRows(end+1,:)={yr,Period{iy},sum(mask), ...
        scores(1),bestNRMSE,models{ii(1)},structPct,hydroPct, ...
        mixedPct,hydroHistory,hydroFinal,overall(1)+overall(2),overall(3)}; %#ok<AGROW>
    fprintf('%d: best=%s (%.4f); structural classes=%.2f%%; hydro=%.2f%%; within hydro historical=%.2f%%\n', ...
        yr,models{ii(1)},scores(ii(1)),structPct,hydroPct,hydroHistory);
end
out=fullfile(cfg.outputRoot,'annual_analysis');
if ~isfolder(out),mkdir(out);end
T=cell2table(modelRows,'VariableNames', ...
    {'TestYear','Period','Model','MeanNRMSE','Rank','CommonCells'});
writetable(T,fullfile(out,'annual_model_rankings.csv'));
T=cell2table(winnerRows,'VariableNames', ...
    {'TestYear','Period','Model','WinnerPct','WinnerCells'});
writetable(T,fullfile(out,'annual_winner_shares.csv'));
T=cell2table(classRows,'VariableNames', ...
    {'TestYear','Period','GainCategory','CategoryPct','CategoryCells'});
writetable(T,fullfile(out,'annual_gain_class_shares.csv'));
T=cell2table(componentRows,'VariableNames', ...
    {'TestYear','Period','GainCategory','Component','ContributionPct', ...
     'WeightedComponentSum','WeightedTotalGainSum'});
writetable(T,fullfile(out,'annual_within_class_components.csv'));
T=cell2table(overallRows,'VariableNames', ...
    {'TestYear','Period','Component','ContributionPct','WeightedComponentSum'});
writetable(T,fullfile(out,'annual_overall_component_shares.csv'));
T=cell2table(summaryRows,'VariableNames', ...
    {'TestYear','Period','CommonCells','ARMeanNRMSE','OracleMeanNRMSE', ...
     'FirstRankedModel','LinearPlusNonlinearClassPct','HydroClassPct', ...
     'MixedClassPct','HistoryShareWithinHydroPct', ...
     'HydroShareWithinHydroPct','OverallStructuralComponentPct', ...
     'OverallHydroComponentPct'});
writetable(T,fullfile(out,'annual_core_summary.csv'));
% DESCRIPTIVE period comparison only; we do not select a six-year winner.
oldMean=mean(cell2mat(summaryRows(4:6,[7:13])),1);
newMean=mean(cell2mat(summaryRows(1:3,[7:13])),1);
statNames={'LinearPlusNonlinearClassPct','HydroClassPct','MixedClassPct', ...
  'HistoryShareWithinHydroPct','HydroShareWithinHydroPct', ...
  'OverallStructuralComponentPct','OverallHydroComponentPct'};
G=table(string(statNames(:)),newMean(:),oldMean(:), ...
    (newMean(:)-oldMean(:)), ...
    'VariableNames',{'Statistic','Mean2019to2021','Mean2022to2024', ...
                     'DifferenceEarlierMinusLater'});
writetable(G,fullfile(out,'period_comparison_descriptive.csv'));
Sample=struct('NRMSE9byYear',E,'testYears',years, ...
   'selectedValidRow',double(P.selectedValidRow(:)), ...
   'selectedLinearIndex',double(P.selectedLinearIndex(:)), ...
   'selectedLat',double(P.selectedLat(:)), ...
   'selectedLon',double(P.selectedLon(:)), ...
   'stratificationWeight',weights,'winner',Winner,'gainClass',Class, ...
   'componentNRMSEGain',Components,'modelNames',{models}, ...
   'linearSourceProvenance',{provenance});
save(fullfile(out,'annual_2000_compact_results.mat'),'-struct','Sample','-v7.3');
SKLGB_ext_period_bootstrap(Sample,out);
fprintf('Analysis output: %s\n',out);
end

function [win,cls,comp,total]=attributeYear(X)
n=size(X,1); win=nan(n,1); cls=nan(n,1); comp=zeros(n,4);
total=zeros(n,1);
for i=1:n
    row=X(i,:);
    if any(~isfinite(row)), continue; end
    [best,w]=min(row);
    [bestL,~]=min(row(1:3));
    [bestHist,ref]=min(row(1:6));
    a=row(1);
    win(i)=w; total(i)=a-best;
    if w==1
        cls(i)=1; % no gain
    elseif w<=3
        cls(i)=2; comp(i,1)=a-best;
    elseif w<=6
        cls(i)=3; comp(i,1)=a-bestL; comp(i,2)=bestL-best;
    elseif ref==w-3
        cls(i)=4; comp(i,1)=a-bestL;
        comp(i,2)=bestL-bestHist; comp(i,3)=bestHist-best;
    elseif ref<=3
        cls(i)=5; comp(i,1)=a-bestHist; comp(i,4)=bestHist-best;
    else
        cls(i)=6; comp(i,1)=a-bestL;
        comp(i,2)=bestL-bestHist; comp(i,4)=bestHist-best;
    end
    if abs(sum(comp(i,:))-total(i))>1e-9
        error('Attribution decomposition failed at row %d',i);
    end
end
end
function s=categoryName(k)
names={'No gain','Linear gain','Nonlinear gain','Hydroclimatic gain','Mixed-L','Mixed-N'};
s=names{k};
end
function s=componentName(k)
names={'Linear','Nonlinear','Hydroclimatic','Mixed'};
s=names{k};
end
