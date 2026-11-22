function SKLGB_ext_period_bootstrap(Sample, out)
% Exploratory paired spatial-block bootstrap for EARLY vs LATE 3-year
% PERIOD-AVERAGED ANNUAL statistics (not a pooled six-year winner).
% These CIs condition on the six observed years: they do NOT quantify
% uncertainty about unobserved years or prove equivalence/stationarity.

w=double(Sample.stratificationWeight(:));
lat=double(Sample.selectedLat(:)); lon=double(Sample.selectedLon(:));
coords=[floor((lat+90)/3),floor((lon+180)/3)];
[~,~,block]=unique(coords,'rows','stable');
nb=max(block); n=length(w); ny=6;
% Weighted block totals for seven core statistics.
BW=zeros(nb,ny); BS=zeros(nb,ny); BH=zeros(nb,ny); BM=zeros(nb,ny);
BHydroTotal=zeros(nb,ny); BHydroHist=zeros(nb,ny); BHydroAdd=zeros(nb,ny);
BAllTotal=zeros(nb,ny); BAllStruct=zeros(nb,ny); BAllHydro=zeros(nb,ny);
for iy=1:ny
    cls=Sample.gainClass(:,iy);
    C=Sample.componentNRMSEGain(:,:,iy);
    complete=isfinite(cls) & all(isfinite(C),2);
    C(~complete,:)=0;
    ww=w; ww(~complete)=0;
    c=cls;
    BW(:,iy)=accumarray(block,ww,[nb,1]);
    BS(:,iy)=accumarray(block,ww.*double(c==2 | c==3),[nb,1]);
    BH(:,iy)=accumarray(block,ww.*double(c==4),[nb,1]);
    BM(:,iy)=accumarray(block,ww.*double(c==5 | c==6),[nb,1]);
    gain=sum(C,2);
    hg=ww.*double(c==4).*gain;
    hh=ww.*double(c==4).*(C(:,1)+C(:,2));
    ha=ww.*double(c==4).*C(:,3);
    BHydroTotal(:,iy)=accumarray(block,hg,[nb,1]);
    BHydroHist(:,iy)=accumarray(block,hh,[nb,1]);
    BHydroAdd(:,iy)=accumarray(block,ha,[nb,1]);
    BAllTotal(:,iy)=accumarray(block,ww.*gain,[nb,1]);
    BAllStruct(:,iy)=accumarray(block,ww.*(C(:,1)+C(:,2)),[nb,1]);
    BAllHydro(:,iy)=accumarray(block,ww.*C(:,3),[nb,1]);
end
nRep=1000;
oldRng=rng; cleanup=onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(20260605,'twister');
Diff=nan(nRep,7);
for rep=1:nRep
    idx=randi(nb,nb,1);
    V=nan(ny,7);
    for iy=1:ny
        den=sum(BW(idx,iy)); hydroG=sum(BHydroTotal(idx,iy));
        gainG=sum(BAllTotal(idx,iy));
        if den<=0 || hydroG<=0 || gainG<=0,continue;end
        V(iy,1)=100*sum(BS(idx,iy))/den;
        V(iy,2)=100*sum(BH(idx,iy))/den;
        V(iy,3)=100*sum(BM(idx,iy))/den;
        V(iy,4)=100*sum(BHydroHist(idx,iy))/hydroG;
        V(iy,5)=100*sum(BHydroAdd(idx,iy))/hydroG;
        V(iy,6)=100*sum(BAllStruct(idx,iy))/gainG;
        V(iy,7)=100*sum(BAllHydro(idx,iy))/gainG;
    end
    Diff(rep,:)=mean(V(1:3,:),1,'omitnan')-mean(V(4:6,:),1,'omitnan');
end
stat={'LinearPlusNonlinearClassPct';'HydroClassPct';'MixedClassPct'; ...
      'HistoryShareWithinHydroPct';'HydroShareWithinHydroPct'; ...
      'OverallStructuralComponentPct';'OverallHydroComponentPct'};
% For MATLAB CSV use readtable instead of load, for compatibility.
ref=readtable(fullfile(out,'period_comparison_descriptive.csv'));
lo=prctile(Diff,2.5,1)'; hi=prctile(Diff,97.5,1)';
T=table(string(stat),ref.DifferenceEarlierMinusLater,lo,hi, ...
    repmat(3,7,1),repmat(nRep,7,1), ...
    'VariableNames',{'Statistic','ObservedDifference','CI95Lower', ...
    'CI95Upper','SpatialBlockDegrees','BootstrapReplicates'});
writetable(T,fullfile(out,'period_comparison_spatial_bootstrap_CI.csv'));
fprintf('Exploratory period-difference block-bootstrap CI: %d resamples.\n',nRep);
end
