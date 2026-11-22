function SKLGB_merge_2000_results()
% Strictly merge the two 1,000-cell outputs of the FIXED 2,000-cell study.
% Copy each computer's compact_sensitivity_results.mat into:
%   merge_inputs/PC1/compact_sensitivity_results.mat
%   merge_inputs/PC2/compact_sensitivity_results.mat
% Then run SKLGB_merge_2000_results on one computer.
% Refuse merge for missing, duplicated, wrongly assigned or incomplete data.

here = fileparts(mfilename('fullpath'));
inputRoot = fullfile(here,'merge_inputs');
outputRoot = fullfile(here,'merged_results');
files = {fullfile(inputRoot,'PC1','compact_sensitivity_results.mat'), ...
         fullfile(inputRoot,'PC2','compact_sensitivity_results.mat')};
for k=1:2
    if ~isfile(files{k})
        error('Missing computer %d compact result: %s',k,files{k});
    end
end
A=load(files{1});
B=load(files{2});
P=load(fullfile(here,'SKLGB_hyperparameter_sample2000.mat'));
if numel(P.selectedValidRow)~=2000 || ...
        ~isequal(double(P.selectedValidRow(:)),sort(double(P.selectedValidRow(:))))
    error('Fixed 2,000-grid sample is invalid.');
end
validateMachine(A,P,1,files{1});
validateMachine(B,P,2,files{2});
if ~isequal(A.models(:),B.models(:)) || ~isequal(A.tagList,B.tagList)
    error('Computer model names or tested hyperparameter tags differ.');
end
if ~isequal(size(A.score),size(B.score)) || ...
        ~isequal(size(A.valScore),size(B.valScore)) || ...
        ~isequal(size(A.bestP),size(B.bestP)) || ...
        ~isequal(size(A.bestX),size(B.bestX))
    error('Computer output dimensions are incompatible.');
end
if ~isequal(size(A.score),[1000 4 7 3])
    error('Expected 1,000-by-4-by-7-by-3 score arrays from both computers.');
end

allRows=[double(A.originalValidRows(:)); double(B.originalValidRows(:))];
[originalValidRows,order]=sort(allRows);
expectedRows=double(P.selectedValidRow(:));
if numel(unique(originalValidRows))~=2000 || ...
        ~isequal(originalValidRows,expectedRows)
    error('Combined grid rows have duplicates, missing cells, or unexpected locations.');
end
originalLinearIndex=concatSort(A,B,'originalLinearIndex',order);
validLat=concatSort(A,B,'validLat',order);
validLon=concatSort(A,B,'validLon',order);
designWeight=concatSort(A,B,'designWeight',order);
stratumId=concatSort(A,B,'stratumId',order);
originalAnnualNRMSE9=concatSort(A,B,'originalAnnualNRMSE9',order);
originalMeanGainClass=concatSort(A,B,'originalMeanGainClass',order);
originalAnnualGainClass=concatSort(A,B,'originalAnnualGainClass',order);
originalAnnualWinnerModel=concatSort(A,B,'originalAnnualWinnerModel',order);
score=concatSort(A,B,'score',order);
valScore=concatSort(A,B,'valScore',order);
bestP=concatSort(A,B,'bestP',order);
bestX=concatSort(A,B,'bestX',order);

if ~isequal(originalLinearIndex,double(P.selectedLinearIndex(:))) || ...
        ~isequal(originalAnnualNRMSE9,double(P.annualNRMSE9)) || ...
        ~isequal(originalMeanGainClass,double(P.meanGainClass(:))) || ...
        ~isequal(originalAnnualGainClass,double(P.annualGainClass)) || ...
        ~isequal(originalAnnualWinnerModel,double(P.annualWinnerModel)) || ...
        ~isequal(stratumId,double(P.stratumId(:))) || ...
        any(abs(validLat-double(P.selectedLat(:)))>1e-10) || ...
        any(abs(validLon-double(P.selectedLon(:)))>1e-10)
    error('Combined cells disagree with the preassigned original 2,000-grid manifest.');
end
pop18=reshape(double(P.populationStratumCounts).',[],1);
n18=accumarray(stratumId,1,[18 1]);
expectedWeight=pop18(stratumId)./n18(stratumId);
if any(abs(expectedWeight-designWeight)>1e-10) || ...
        abs(sum(designWeight)-17284)>1e-7
    error('Merged design weights are inconsistent with the full 18 strata.');
end
if any(~isfinite(score(:))) || any(~isfinite(valScore(:, :, 2:7, :)),'all') || ...
        any(~isfinite(bestP(:,:,2:7,:)),'all')
    error('Merged data contain missing or nonfinite completed model metrics.');
end
models=A.models;
tagList=A.tagList;
wasFit=A.wasFit & B.wasFit;
computeMinutesPerPC=cat(3,double(A.computeMinutes),double(B.computeMinutes));
successCountPerPC=cat(3,double(A.successCount),double(B.successCount));
sourceComputer=[ones(1000,1);2*ones(1000,1)];
sourceComputer=sourceComputer(order);
if ~isequal(sourceComputer(1:2:end),ones(1000,1)) || ...
        ~isequal(sourceComputer(2:2:end),2*ones(1000,1))
    error('Merged machine allocation differs from the frozen alternating partition.');
end
if ~isfolder(outputRoot),mkdir(outputRoot);end
out=fullfile(outputRoot,'SKLGB_merged_2000_results.mat');
save(out,'originalValidRows','originalLinearIndex','validLat','validLon', ...
    'designWeight','stratumId','originalAnnualNRMSE9', ...
    'originalMeanGainClass','originalAnnualGainClass', ...
    'originalAnnualWinnerModel','score','valScore','bestP','bestX', ...
    'wasFit','computeMinutesPerPC','successCountPerPC', ...
    'models','tagList','sourceComputer','-v7.3');

% Compact computation-cost records from the actually completed programs.
rows=cell(24,7);
k=0;
for m=1:4
    for j=2:7
        k=k+1;
        c1=computeMinutesPerPC(m,j,1);
        c2=computeMinutesPerPC(m,j,2);
        rows(k,:)={models{m},tagList{m,j},c1,c2,c1+c2, ...
            (c1+c2)*(17284/2000)/60, ...
            max(c1,c2)/60};
    end
end
R=cell2table(rows,'VariableNames',{'Model','Configuration', ...
    'PC1ComputeMinutes','PC2ComputeMinutes','TotalTwoPCComputeMinutes', ...
    'ExtrapolatedGlobalSinglePC_Hours','TwoPCStageLongest_Hours'});
writetable(R,fullfile(outputRoot,'SKLGB_runtime_2000.csv'));
T=table(originalValidRows,originalLinearIndex,validLat,validLon, ...
    stratumId,designWeight,sourceComputer,'VariableNames', ...
    {'OriginalValidRow','LinearIndex','Latitude','Longitude', ...
    'Stratum','DesignWeight','MachineIndex'});
writetable(T,fullfile(outputRoot,'SKLGB_merged_2000_manifest.csv'));
fprintf('\nSTRICT MERGE PASSED: 2,000 unique cells = PC1 1,000 + PC2 1,000.\n');
fprintf('All 24 perturbation configurations have 2,000 x 3 successful tests.\n');
fprintf('Merged MAT: %s\n',out);
SKLGB_analyze_merged_2000(out,outputRoot);
end

function validateMachine(S,P,computer,path)
required={'originalValidRows','originalLinearIndex','validLat','validLon', ...
    'designWeight','stratumId','originalAnnualNRMSE9', ...
    'originalMeanGainClass','originalAnnualGainClass', ...
    'originalAnnualWinnerModel','models','tagList','score','valScore', ...
    'bestP','bestX','wasFit','computeMinutes','successCount', ...
    'expectedGridCount','expectedMachineCount'};
for k=1:numel(required)
    if ~isfield(S,required{k})
        error('PC%d result is missing variable %s: %s',computer,required{k},path);
    end
end
if S.expectedGridCount~=2000 || S.expectedMachineCount~=2
    error('PC%d output was not created by fixed 2,000-grid two-PC design.',computer);
end
expected=double(P.selectedValidRow(computer:2:end));
if ~isequal(double(S.originalValidRows(:)),expected(:)) || ...
        numel(S.originalValidRows)~=1000
    error('PC%d did not process its exact preallocated 1,000 grid cells.',computer);
end
if ~isequal(size(S.wasFit),[4,7]) || ~all(all(S.wasFit(:,2:7))) || ...
        ~isequal(size(S.successCount),[4,7]) || ...
        any(any(double(S.successCount(:,2:7))~=3000)) || ...
        any(~isfinite(S.computeMinutes(:,2:7)),'all')
    error('PC%d has incomplete, failed, or missing perturbation configurations.',computer);
end
expectedNames={'RF';'RFX';'LSTM';'LSTMX'};
if ~isequal(S.models(:),expectedNames)
    error('PC%d has incorrect model order.',computer);
end
expectedRF={'baseline','trees_025','trees_100','trees_200', ...
    'leaf_001','leaf_010','leaf_020'};
expectedLSTM={'baseline','hidden_016','hidden_064','hidden_128', ...
    'lr_0p001','lr_0p003','lr_0p030'};
if ~isequal(S.tagList(1,:),expectedRF) || ...
        ~isequal(S.tagList(2,:),expectedRF) || ...
        ~isequal(S.tagList(3,:),expectedLSTM) || ...
        ~isequal(S.tagList(4,:),expectedLSTM)
    error('PC%d has incorrect perturbation tag/order.',computer);
end
if any(~isfinite(S.score(:))) || any(~isfinite(S.valScore(:,:,2:7,:)),'all') || ...
        any(~isfinite(S.bestP(:,:,2:7,:)),'all')
    error('PC%d has missing/nonfinite scores or selected lags.',computer);
end
if any(~ismember(S.bestP(:,:,2:7,:),[1,6,12,24,36,48,60]),'all')
    error('PC%d contains a selected lag outside the prescribed p search.',computer);
end
if any(~ismember(S.bestX(:,2,2:7,:),[1,3,6]),'all')
    error('PC%d contains a selected RFX exogenous lag outside x search.',computer);
end
end

function Y=concatSort(A,B,field,order)
X=cat(1,A.(field),B.(field));
ix=repmat({':'},1,ndims(X));
ix{1}=order;
Y=X(ix{:});
Y=double(Y);
end
