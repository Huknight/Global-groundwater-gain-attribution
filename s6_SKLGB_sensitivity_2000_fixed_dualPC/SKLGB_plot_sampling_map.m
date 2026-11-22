function SKLGB_plot_sampling_map(cfg)
% Plot selected sample distribution in the same plain global-map style as
% s1_overall_performance3.m: 1-degree grid, Times New Roman, 300-dpi PNG.
% Fig_sample_global: six attribution colors for all requested cells.
% Fig_sample_machine: allocation to current computer (two-PC check).
% Uses original full land/data mask, not an approximate coastline.

P = load(cfg.sampleFile, 'selectedValidRow','selectedLinearIndex', ...
    'selectedLon','selectedLat','meanGainClass','stratumId');
if any(strcmp('nestedPriority',who('-file',cfg.sampleFile)))
    T = load(cfg.sampleFile,'nestedPriority');
    positions = find(double(T.nestedPriority(:)) <= cfg.nPilot);
else
    positions = unique(round(linspace(1,numel(P.selectedValidRow),cfg.nPilot)))';
end
if numel(positions) ~= cfg.nPilot
    error('Sample map: invalid total sample count.');
end
positions = positions(:);
localAssignment = mod((0:numel(positions)-1)',cfg.machineCount)+1;
isLocal = localAssignment == cfg.machineIndex;

D = load(cfg.originalMasterFile,'CommonValid','lat','lon');
validMask = logical(D.CommonValid);
lat = double(D.lat(:));
lon = double(D.lon(:));
if ~isequal(size(validMask),[150,360]) || nnz(validMask)~=17284
    error('Sample map: unexpected common valid mask.');
end
[Lon2D,Lat2D] = meshgrid(lon,lat);
ix = double(P.selectedLinearIndex(positions));
if any(~validMask(ix)) || numel(unique(ix))~=numel(ix)
    error('Sample map: sample includes invalid or duplicate grid cells.');
end
sLat = double(P.selectedLat(positions));
sLon = double(P.selectedLon(positions));
sClass = double(P.meanGainClass(positions));
if any(abs(Lat2D(ix)-sLat)>1e-8) || any(abs(Lon2D(ix)-sLon)>1e-8)
    error('Sample map: sample-coordinate mismatch.');
end
classNames = {'No gain','Linear','Nonlinear','Hydroclimatic','Mixed-L','Mixed-N'};
classColors = [0.52 0.52 0.52; 0.84 0.14 0.14; 0.10 0.38 0.88; ...
    0.07 0.62 0.27; 0.99 0.69 0.05; 0.62 0.30 0.78];

% Fig. 1: six original three-year mean gain-attribution classes
fig = figure('Color','w','Position',[100 100 1200 650], ...
    'Visible','off');
ax = axes('Parent',fig);
plotBase(ax,Lon2D,Lat2D,validMask);
hold(ax,'on');
leg = gobjects(0); legTxt = {};
for c=1:6
    subset = sClass==c;
    if nnz(subset)==0,continue;end
    h = scatter(ax,sLon(subset),sLat(subset), ...
        chooseSize(cfg.nPilot),classColors(c,:),'filled', ...
        'MarkerEdgeColor',[0.18 0.18 0.18], ...
        'LineWidth',0.2);
    leg(end+1)=h; %#ok<AGROW>
    legTxt{end+1}=sprintf('%s (n=%d)',classNames{c},nnz(subset)); %#ok<AGROW>
end
formatAxes(ax,sprintf('Hyperparameter sensitivity sample: %d of 17,284 grid cells',cfg.nPilot));
if ~isempty(leg)
    legend(ax,leg,legTxt,'Location','southoutside', ...
        'Orientation','horizontal','NumColumns',3,'FontName','Times New Roman', ...
        'FontSize',10,'Box','off');
end
exportgraphics(fig,fullfile(cfg.outputRoot,'Fig_sample_global.png'),'Resolution',300);
close(fig);

% Fig. 2: which points this computer will calculate
fig = figure('Color','w','Position',[100 100 1200 650], ...
    'Visible','off');
ax = axes('Parent',fig);
plotBase(ax,Lon2D,Lat2D,validMask);
hold(ax,'on');
other = scatter(ax,sLon(~isLocal),sLat(~isLocal), ...
    chooseSize(cfg.nPilot),[0.90 0.61 0.20],'filled', ...
    'MarkerEdgeColor','none');
current = scatter(ax,sLon(isLocal),sLat(isLocal), ...
    chooseSize(cfg.nPilot),[0.12 0.38 0.80],'filled', ...
    'MarkerEdgeColor','none');
formatAxes(ax,sprintf('Computer %d of %d: %d selected cells of %d', ...
    cfg.machineIndex,cfg.machineCount,nnz(isLocal),cfg.nPilot));
legend(ax,[current,other],{sprintf('This computer (n=%d)',nnz(isLocal)), ...
    sprintf('Other computer(s) (n=%d)',nnz(~isLocal))}, ...
    'Location','southoutside','Orientation','horizontal', ...
    'FontName','Times New Roman','FontSize',11,'Box','off');
exportgraphics(fig,fullfile(cfg.outputRoot,'Fig_sample_machine.png'),'Resolution',300);
close(fig);

% Machine allocation audit table, also useful before starting long jobs.
T=table((1:numel(positions))',double(P.selectedValidRow(positions)), ...
    ix,sLat,sLon,sClass, ...
    double(P.stratumId(positions)),localAssignment, ...
    'VariableNames',{'SampleRow','OriginalValidRow','LinearIndex', ...
    'Latitude','Longitude','MeanGainClass','Stratum','MachineIndex'});
writetable(T,fullfile(cfg.outputRoot,'sample_machine_allocation.csv'));
fprintf('Sample maps saved to %s (n=%d; local=%d).\n', ...
    cfg.outputRoot,cfg.nPilot,nnz(isLocal));
end

function plotBase(ax,Lon2D,Lat2D,validMask)
base = double(validMask);
base(~validMask) = NaN;
pcolor(ax,Lon2D,Lat2D,base);
shading(ax,'flat');
colormap(ax,[0.80 0.80 0.80]);
clim(ax,[0,1]);
set(ax,'Color',[1 1 1]);
end

function formatAxes(ax,figureTitle)
set(ax,'YDir','normal','FontName','Times New Roman','FontSize',12, ...
    'LineWidth',0.8,'Box','on');
axis(ax,'equal');
xlim(ax,[-180,180]);
ylim(ax,[-60,90]);
xticks(ax,-180:60:180);
yticks(ax,-60:30:90);
xlabel(ax,'Longitude','FontName','Times New Roman','FontSize',13);
ylabel(ax,'Latitude','FontName','Times New Roman','FontSize',13);
title(ax,figureTitle,'FontName','Times New Roman', ...
    'FontSize',16,'FontWeight','normal');
end

function s=chooseSize(n)
if n<=10, s=42;
elseif n<=400,s=22;
else,s=11;
end
end
