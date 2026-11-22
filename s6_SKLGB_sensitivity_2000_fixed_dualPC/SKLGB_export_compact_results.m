function SKLGB_export_compact_results(outputRoot)
% Export all finished scenarios as one small file for transfer and merging.
% It is not necessary to upload large forecast MAT files for first analysis.
% Call after completion, or manually following an interrupted run.
M = load(fullfile(outputRoot,'pilot10_master.mat'), ...
    'pilotBaselineNRMSE9','pilotMeanGainClass','pilotAnnualGainClass', ...
    'pilotAnnualWinnerModel','pilotDesignWeight','pilotStratumId', ...
    'originalValidRows','originalLinearIndex','validLat','validLon');
n = numel(M.originalValidRows);
models = {'RF','RFX','LSTM','LSTMX'};
rfTags={'baseline','trees_025','trees_100','trees_200', ...
    'leaf_001','leaf_010','leaf_020'};
lstmTags={'baseline','hidden_016','hidden_064','hidden_128', ...
    'lr_0p001','lr_0p003','lr_0p030'};
modelFiles={'RF_univariate_global_200401_202412.mat', ...
    'RFX_multivariate_global_200401_202412.mat', ...
    'LSTM_univariate_global_200401_202412.mat', ...
    'LSTMX_multivariate_global_200401_202412.mat'};

score=nan(n,4,7,3);
valScore=nan(n,4,7,3);
bestP=nan(n,4,7,3);
bestX=nan(n,4,7,3);
wasFit=false(4,7);
computeMinutes=nan(4,7);
successCount=zeros(4,7);
tagList=cell(4,7);
sourceModelColumn=[5,8,6,9];
for m=1:4
    model=models{m};
    if m<=2,ts=rfTags;else,ts=lstmTags;end
    for j=1:7
        tag=ts{j};
        tagList{m,j}=tag;
        if j==1
            score(:,m,j,:)=reshape(M.pilotBaselineNRMSE9(:,sourceModelColumn(m),:), ...
                [n,1,1,3]);
            successCount(m,j)=3*n;
        end
        f=fullfile(outputRoot,model,tag,modelFiles{m});
        if ~isfile(f), continue; end
        vars={[model '_NRMSE_valid'],[model '_Val_NRMSE_best_valid'], ...
            [model '_best_p_valid'],[model '_TestSuccess_valid'],'Meta'};
        if strcmp(model,'RFX'),vars{end+1}='RFX_best_x_valid';end
        D=load(f,vars{:});
        if ~isequal(double(D.Meta.pilotOriginalValidRows(:)), ...
                double(M.originalValidRows(:)))
            error('Compact exporter: grid identity mismatch for %s %s.',model,tag);
        end
        tmp=double(D.([model '_NRMSE_valid']));
        if ~isequal(size(tmp),[n,3]),error('Bad NRMSE shape: %s',f);end
        score(:,m,j,:)=reshape(tmp,[n,1,1,3]);
        valScore(:,m,j,:)=reshape(double(D.([model '_Val_NRMSE_best_valid'])), ...
            [n,1,1,3]);
        bestP(:,m,j,:)=reshape(double(D.([model '_best_p_valid'])), ...
            [n,1,1,3]);
        if strcmp(model,'RFX')
            bestX(:,m,j,:)=reshape(double(D.RFX_best_x_valid),[n,1,1,3]);
        end
        successCount(m,j)=nnz(D.([model '_TestSuccess_valid']));
        computeMinutes(m,j)=double(D.Meta.computeTimeMinutes);
        wasFit(m,j)=true;
    end
end
out=fullfile(outputRoot,'compact_sensitivity_results.mat');
% The master is machine-specific. The fixed sample manifest is checked
% again when both computer outputs are merged.
% Store the assignment by reference to the actual local original rows.
expectedGridCount=2000;
expectedMachineCount=2;
originalValidRows=double(M.originalValidRows(:));
originalLinearIndex=double(M.originalLinearIndex(:));
validLat=double(M.validLat(:));
validLon=double(M.validLon(:));
designWeight=double(M.pilotDesignWeight(:));
stratumId=double(M.pilotStratumId(:));
originalAnnualNRMSE9=double(M.pilotBaselineNRMSE9);
originalMeanGainClass=double(M.pilotMeanGainClass(:));
originalAnnualGainClass=double(M.pilotAnnualGainClass);
originalAnnualWinnerModel=double(M.pilotAnnualWinnerModel);
save(out,'originalValidRows','originalLinearIndex','validLat','validLon', ...
    'designWeight','stratumId','originalAnnualNRMSE9', ...
    'originalMeanGainClass','originalAnnualGainClass', ...
    'originalAnnualWinnerModel','models','tagList','score','valScore', ...
    'bestP','bestX','wasFit','computeMinutes','successCount', ...
    'expectedGridCount','expectedMachineCount','-v7.3');
fprintf('Compact result saved: %s (%d cells, %d trained config files).\n', ...
    out,n,nnz(wasFit));
end
