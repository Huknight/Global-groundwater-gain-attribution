function SKLGB_mean_analyze(cfg)
% Paired comparisons on EXACTLY the original 2000 cells; no global extrapolation.
% 2022/23/24 annual + PRIMARY 2022-24 mean-NRMSE attribution.
if nargin<1,cfg=SKLGB_mean_config();end
P=load(cfg.sampleFile,'selectedValidRow','selectedLinearIndex', ...
    'selectedLat','selectedLon','designWeight','annualNRMSE9', ...
    'annualGainClass','annualWinnerModel','meanGainClass');
models={'AR','ARIMA','SARIMA','SVM','RF','LSTM','SVMX','RFX','LSTMX'};
classes={'No gain','Linear gain','Nonlinear gain','Hydroclimatic gain','Mixed-L','Mixed-N'};
components={'Linear','Nonlinear','Hydroclimatic','Mixed'};
n=2000;W=double(P.designWeight(:));
assert(isequal(size(P.annualNRMSE9),[n,9,3]) && ...
    all(W>0) && abs(sum(W)-17284)<1e-6);
ref=double(P.annualNRMSE9);
alt=nan(n,9,3);
for im=1:9
    name=models{im};
    f=fullfile(cfg.outputRoot,name,sprintf('%s_monthmean_2022_2024.mat',name));
    if ~isfile(f),error('Missing %s model output: %s',name,f);end
    S=load(f,[name '_NRMSE_valid'],[name '_TestSuccess_valid'], ...
        'validLinearIndex','testYears','Meta');
    Z=double(S.([name '_NRMSE_valid']));
    success=S.([name '_TestSuccess_valid']);
    assert(isequal(size(Z),[n,3]) && isequal(size(success),[n,3]));
    assert(isequal(double(S.validLinearIndex(:)), ...
        double(P.selectedLinearIndex(:))) && ...
        isequal(double(S.testYears(:))',[2022 2023 2024]) && ...
        S.Meta.randomSeed==20260605 && ...
        strcmp(S.Meta.targetDefinition, ...
          'difference of consecutive 1-degree monthly-mean GWS states'));
    Z(~success)=NaN;
    alt(:,im,:)=reshape(Z,n,1,3);
end
% Verify the original attribution implementation matches stored 2000 labels.
for iy=1:3
    [w,c]=attribute(ref(:,:,iy));
    assert(all(w==double(P.annualWinnerModel(:,iy))) && ...
        all(c==double(P.annualGainClass(:,iy))), ...
        'Original gain-category / winner reproduction check failed.');
end
[~,baseMeanCat]=attribute(mean(ref,3));
assert(all(baseMeanCat==double(P.meanGainClass(:))), ...
    'Primary original 3-year mean gain categories cannot be reproduced.');

out=fullfile(cfg.outputRoot,'analysis_monthmean_vs_monthend');
if ~isfolder(out),mkdir(out);end
RankRows=cell(0,7);WinRows=cell(0,7);CatRows=cell(0,7);
CompRows=cell(0,8);OverallRows=cell(0,7);
SummaryRows=cell(0,12);Pairs=cell(0,9);Switch=cell(0,7);
Result=struct();Result.OriginalNRMSE=ref;Result.MonthMeanNRMSE=alt;
Result.modelNames=models;Result.gridRows=double(P.selectedValidRow(:));
Result.linearIndex=double(P.selectedLinearIndex(:));
Result.latitude=double(P.selectedLat(:));Result.longitude=double(P.selectedLon(:));
Result.stratificationWeights=W;
Result.winner=zeros(n,2,4);Result.gainClass=zeros(n,2,4);
Result.gainComponents=nan(n,4,2,4);
for iy=1:4
    if iy<=3
        period=sprintf('%d',2021+iy);
        origX=ref(:,:,iy); altX=alt(:,:,iy);
    else
        period='2022-2024 mean';
        origX=mean(ref,3);altX=mean(alt,3);% strict common mask computed below.
    end
    if iy<=3
        valid=all(isfinite(origX),2)&all(isfinite(altX),2);
    else
        valid=squeeze(all(all(isfinite(ref)&isfinite(alt),2),3));
    end
    if nnz(valid)<0.95*n
        error('%s: only %d/%d complete paired nine-model grid cells. Inspect failed models.', ...
            period,nnz(valid),n);
    end
    if nnz(valid)<n
        warning('SKLGB:IncompleteSample', ...
            '%s: %d complete cases out of 2000; all comparisons use matched cells only.', ...
            period,nnz(valid));
    end
    pairedW=W(valid);totalW=sum(pairedW);
    winOriginal=nan(n,1);winAlternative=nan(n,1);
    clsOriginal=nan(n,1);clsAlternative=nan(n,1);
    for target=1:2
        X=origX;
        targetName='Month-end dGWS (original)';
        if target==2
            X=altX;targetName='Monthly-mean dGWS (alternative)';
        end
        [win,cls,comp,gain]=attribute(X);
        win(~valid)=NaN;cls(~valid)=NaN;comp(~valid,:)=NaN;
        Result.winner(:,target,iy)=win;
        Result.gainClass(:,target,iy)=cls;
        Result.gainComponents(:,:,target,iy)=comp;
        if target==1,winOriginal=win;clsOriginal=cls;
        else,winAlternative=win;clsAlternative=cls;end
        means=nan(1,9);
        for m=1:9
            means(m)=sum(W(valid).*X(valid,m))/totalW;
        end
        [~,order]=sort(means);
        for m=1:9
            RankRows(end+1,:)={period,targetName,models{m}, ...
                means(m),find(order==m,1),nnz(valid),totalW}; %#ok<AGROW>
            WinRows(end+1,:)={period,targetName,models{m}, ...
                100*sum(W(valid&win==m))/totalW,nnz(valid&win==m), ...
                nnz(valid),totalW}; %#ok<AGROW>
        end
        for c=1:6
            mask=valid&cls==c;
            p=100*sum(W(mask))/totalW;
            CatRows(end+1,:)={period,targetName,classes{c},p,nnz(mask), ...
                nnz(valid),totalW}; %#ok<AGROW>
            if c==1,continue;end
            tot=sum(W(mask).*gain(mask));
            for k=1:4
                val=sum(W(mask).*comp(mask,k));
                share=NaN;if tot>0,share=100*val/tot;end
                CompRows(end+1,:)={period,targetName,classes{c}, ...
                    components{k},share,val,tot,nnz(mask)}; %#ok<AGROW>
            end
        end
        gainTot=sum(W(valid).*gain(valid));
        for k=1:4
            share=NaN;if gainTot>0,share=100*sum(W(valid).*comp(valid,k))/gainTot;end
            OverallRows(end+1,:)={period,targetName,components{k}, ...
                share,sum(W(valid).*comp(valid,k)),gainTot,nnz(valid)}; %#ok<AGROW>
        end
        hydro=valid&cls==4;
        hTot=sum(W(hydro).*gain(hydro));
        hist=NaN;finalH=NaN;
        if hTot>0
            hist=100*sum(W(hydro).*(comp(hydro,1)+comp(hydro,2)))/hTot;
            finalH=100*sum(W(hydro).*comp(hydro,3))/hTot;
        end
        structP=100*sum(W(valid&(cls==2|cls==3)))/totalW;
        hydroP=100*sum(W(hydro))/totalW;
        mixedP=100*sum(W(valid&(cls==5|cls==6)))/totalW;
        SummaryRows(end+1,:)={period,targetName,nnz(valid),totalW, ...
            models{order(1)},means(order(1)), ...
            100*sum(W(valid&win==8))/totalW, ...
            structP,hydroP,mixedP,hist,finalH}; %#ok<AGROW>
    end
    winnerAgree=100*sum(W(valid&(winOriginal==winAlternative)))/totalW;
    classAgree=100*sum(W(valid&(clsOriginal==clsAlternative)))/totalW;
    Pairs(end+1,:)={period,nnz(valid),totalW,winnerAgree,classAgree, ...
        100-winnerAgree,100-classAgree, ...
        nnz(valid&(clsOriginal==4)&(clsAlternative==4)), ...
        nnz(valid&(clsOriginal==4 | clsAlternative==4))}; %#ok<AGROW>
    for i=1:6
        for j=1:6
            keep=valid&clsOriginal==i&clsAlternative==j;
            Switch(end+1,:)={period,classes{i},classes{j}, ...
                100*sum(W(keep))/totalW,nnz(keep),nnz(valid),totalW}; %#ok<AGROW>
        end
    end
    fprintf('%s: %d paired cells; winner agreement %.2f%%; class agreement %.2f%%.\n', ...
        period,nnz(valid),winnerAgree,classAgree);
end
writetable(cell2table(RankRows,'VariableNames', ...
    {'Period','TargetDefinition','Model','MeanNRMSE','ModelRank','CommonCells','SumDesignWeights'}), ...
    fullfile(out,'model_rankings_per_target.csv'));
writetable(cell2table(WinRows,'VariableNames', ...
    {'Period','TargetDefinition','Model','WinnerPct','WinnerCells','CommonCells','SumDesignWeights'}), ...
    fullfile(out,'winner_proportions_per_target.csv'));
writetable(cell2table(CatRows,'VariableNames', ...
    {'Period','TargetDefinition','GainCategory','CategoryPct','CategoryCells','CommonCells','SumDesignWeights'}), ...
    fullfile(out,'gain_category_proportions_per_target.csv'));
writetable(cell2table(CompRows,'VariableNames', ...
    {'Period','TargetDefinition','GainCategory','Component','ContributionPct', ...
     'WeightedComponentSum','WeightedTotalGainSum','CategoryCells'}), ...
    fullfile(out,'within_category_gain_contributions.csv'));
writetable(cell2table(OverallRows,'VariableNames', ...
    {'Period','TargetDefinition','Component','ContributionPct','WeightedComponentSum', ...
     'WeightedTotalGainSum','CommonCells'}), ...
    fullfile(out,'overall_gain_contributions.csv'));
writetable(cell2table(SummaryRows,'VariableNames', ...
    {'Period','TargetDefinition','CommonCells','SumDesignWeights', ...
     'BestModel','BestModelMeanNRMSE','RFXWinnerPct', ...
     'LinearPlusNonlinearGainPct','HydroclimaticGainPct','MixedGainPct', ...
     'HistoryComponentsWithinHydroPct','HydroComponentWithinHydroPct'}), ...
    fullfile(out,'core_summary_original_vs_mean.csv'));
writetable(cell2table(Pairs,'VariableNames', ...
    {'Period','CommonCells','SumDesignWeights','WinnerAgreementPct', ...
     'GainCategoryAgreementPct','WinnerChangePct','GainCategoryChangePct', ...
     'HydroToHydroCells','HydroInEitherTargetCells'}), ...
    fullfile(out,'paired_gridcell_agreement.csv'));
writetable(cell2table(Switch,'VariableNames', ...
    {'Period','OriginalGainCategory','AlternativeGainCategory','WeightedPctOfSample', ...
     'CellCount','CommonCells','SumDesignWeights'}), ...
    fullfile(out,'gain_category_transition_matrix.csv'));
Result.targetNames={'Original month-end','Alternative monthly mean'};
Result.periodNames={'2022','2023','2024','2022-2024 mean'};
Result.NRMSECrossTargetNote='Absolute NRMSE scales are target-dependent; compare ranking and attribution, not absolute values across targets.';
Result.inference='2000-cell weighted descriptive paired comparison only. No full-domain extrapolation.';
save(fullfile(out,'paired_2000_monthmean_sensitivity.mat'),'-struct','Result','-v7.3');
fprintf('All comparison tables and source arrays saved: %s\n',out);
end

function [win,cls,comp,total]=attribute(X)
n=size(X,1);win=nan(n,1);cls=nan(n,1);
comp=zeros(n,4);total=zeros(n,1);
for i=1:n
    z=X(i,:);if any(~isfinite(z)),continue;end
    [v,w]=min(z);L=min(z(1:3));[H,r]=min(z(1:6));
    ar=z(1);win(i)=w;total(i)=ar-v;
    if w==1
        cls(i)=1;
    elseif w<=3
        cls(i)=2;comp(i,1)=ar-v;
    elseif w<=6
        cls(i)=3;comp(i,1)=ar-L;comp(i,2)=L-v;
    elseif r==w-3
        cls(i)=4;comp(i,1)=ar-L;comp(i,2)=L-H;comp(i,3)=H-v;
    elseif r<=3
        cls(i)=5;comp(i,1)=ar-H;comp(i,4)=H-v;
    else
        cls(i)=6;comp(i,1)=ar-L;comp(i,2)=L-H;comp(i,4)=H-v;
    end
    if abs(sum(comp(i,:))-total(i))>1e-9
        error('Attribution decomposition mismatch at row %d.',i);
    end
end
end
