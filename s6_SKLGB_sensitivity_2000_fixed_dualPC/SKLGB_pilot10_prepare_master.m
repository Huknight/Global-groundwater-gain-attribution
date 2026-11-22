function SKLGB_pilot10_prepare_master(cfg)
% Create a cell-subset master database without changing the spatial grid.
% The original 17,284-row identities are stored for deterministic fit seeds.
% All six predictor time series and all original months are preserved.
if ~isfile(cfg.sampleFile)
    error('Sample file not found: %s', cfg.sampleFile);
end
if ~isfile(cfg.originalMasterFile)
    error('Full master database not found: %s', cfg.originalMasterFile);
end
P = load(cfg.sampleFile);  % Includes nestedPriority in 2000-grid design.
nAvailable = numel(P.selectedValidRow);
if nAvailable~=2000 || cfg.nPilot~=2000 || cfg.machineCount~=2 || ...
        ~isfield(P,'nestedPriority')
    error('The formal version requires the fixed 2,000-grid sample and two PCs.');
end
if cfg.nPilot > nAvailable || cfg.nPilot < 1 || ...
        cfg.nPilot ~= floor(cfg.nPilot)
    error('nPilot must be between 1 and %d.', nAvailable);
end
% First choose the total requested sample. sample2000 nestedPriority
% guarantees n=10 subset n=400 subset n=2000, retaining previous pilot cells.
if isfield(P,'nestedPriority')
    totalPositions = find(double(P.nestedPriority(:)) <= cfg.nPilot);
else
    totalPositions = unique(round(linspace(1,nAvailable,cfg.nPilot)))';
end
if numel(totalPositions) ~= cfg.nPilot
    error('Total sample selection is invalid.');
end
% Split across computers using alternating positions in the full sorted
% sample: no duplicate grid and almost equal local counts.
localAssignment = mod((0:numel(totalPositions)-1)',cfg.machineCount)+1;
selectedSamplePositions = totalPositions(localAssignment == cfg.machineIndex);
if numel(selectedSamplePositions) ~= cfg.nLocal
    error('Per-computer sample size mismatch: expected %d, got %d.', ...
        cfg.nLocal,numel(selectedSamplePositions));
end
originalValidRows = double(P.selectedValidRow(selectedSamplePositions));
originalLinearIndex = double(P.selectedLinearIndex(selectedSamplePositions));
expectedLat = double(P.selectedLat(selectedSamplePositions));
expectedLon = double(P.selectedLon(selectedSamplePositions));
pilotBaselineNRMSE9 = double(P.annualNRMSE9(selectedSamplePositions,:,:));
pilotMeanGainClass = double(P.meanGainClass(selectedSamplePositions));
pilotAnnualGainClass = double(P.annualGainClass(selectedSamplePositions,:));
pilotAnnualWinnerModel = double(P.annualWinnerModel(selectedSamplePositions,:));
% Recompute stratum expansion weights for the *requested* total sample;
% the precomputed file weights correspond to all 2000 cells and must not
% silently be reused for nPilot=10/400.
pop18 = reshape(double(P.populationStratumCounts).',[],1);
allStrata = double(P.stratumId(totalPositions));
count18 = accumarray(allStrata,1,[18,1]);
localStrata = double(P.stratumId(selectedSamplePositions));
pilotDesignWeight = pop18(localStrata) ./ count18(localStrata);
pilotStratumId = double(P.stratumId(selectedSamplePositions));

if any(diff(originalValidRows) <= 0) || any(diff(originalLinearIndex) <= 0)
    error('Selected original row indices must be strictly increasing.');
end
if isfile(cfg.pilotMasterFile) && ~cfg.rebuildPilotMaster
    Q = load(cfg.pilotMasterFile, 'originalValidRows','originalLinearIndex', ...
        'pilotBaselineNRMSE9');
    if isequal(Q.originalValidRows, originalValidRows) && ...
            isequal(Q.originalLinearIndex, originalLinearIndex) && ...
            isequaln(Q.pilotBaselineNRMSE9, pilotBaselineNRMSE9)
        fprintf('Reusing validated pilot master: %s\n', cfg.pilotMasterFile);
        return;
    else
        error(['Existing pilot master does not match the present selection. ', ...
            'Set rebuildPilotMaster=true after clearing old pilot outputs.']);
    end
end

fprintf('Loading full GLDAS master to select %d local cells (of %d total)...\n', ...
    cfg.nLocal,cfg.nPilot);
vars = {'dGWS_valid', 'Prec_valid','Evap_valid','Runoff_valid', ...
    'dSM2m_valid','SurfT_valid','CommonValid','validLinearIndex', ...
    'validLat','validLon','lat','lon','ym'};
D = load(cfg.originalMasterFile, vars{:});
if ~isequal(double(D.validLinearIndex(originalValidRows(:))), originalLinearIndex(:))
    error('Sample rows do not match the full database ordering.');
end
if any(abs(double(D.validLat(originalValidRows(:))) - expectedLat(:)) > 1e-8) || ...
        any(abs(double(D.validLon(originalValidRows(:))) - expectedLon(:)) > 1e-8)
    error('Sample coordinates do not match the full database.');
end
if numel(D.validLinearIndex) ~= 17284 || ...
        ~isequal(size(D.CommonValid), [150, 360])
    error('Unexpected full master mask dimensions / valid-cell count.');
end

% The compact master keeps the original 150-by-360 maps and 1-based
% spatial linear indices. Only the six per-cell time-series matrices shrink.
dGWS_valid = double(D.dGWS_valid(originalValidRows,:));
Prec_valid = double(D.Prec_valid(originalValidRows,:));
Evap_valid = double(D.Evap_valid(originalValidRows,:));
Runoff_valid = double(D.Runoff_valid(originalValidRows,:));
dSM2m_valid = double(D.dSM2m_valid(originalValidRows,:));
SurfT_valid = double(D.SurfT_valid(originalValidRows,:));
lat = double(D.lat(:));
lon = double(D.lon(:));
ym = double(D.ym(:));
validLinearIndex = originalLinearIndex(:);
validLat = expectedLat(:);
validLon = expectedLon(:);
CommonValid = false(150,360);
CommonValid(validLinearIndex) = true;
CommonMask = nan(150,360);
CommonMask(validLinearIndex) = 1;
clear D
inputBlocks = {dGWS_valid,Prec_valid,Evap_valid,Runoff_valid,dSM2m_valid,SurfT_valid};
for k = 1:numel(inputBlocks)
    if ~isequal(size(inputBlocks{k}), [cfg.nLocal,252]) || ...
            any(~isfinite(inputBlocks{k}(:)))
        error('Input block %d has incorrect dimensions or nonfinite values.', k);
    end
end
if ~isequal(validLinearIndex, find(CommonValid))
    error('Pilot spatial index ordering does not match the mask.');
end
[Lon2D,Lat2D] = meshgrid(lon,lat);
if any(abs(Lat2D(validLinearIndex)-validLat)>1e-8) || ...
        any(abs(Lon2D(validLinearIndex)-validLon)>1e-8)
    error('Pilot coordinate mapping failed.');
end

save(cfg.pilotMasterFile, 'dGWS_valid','Prec_valid','Evap_valid', ...
    'Runoff_valid','dSM2m_valid','SurfT_valid', ...
    'CommonValid','CommonMask','validLinearIndex','validLat','validLon', ...
    'lat','lon','ym','originalValidRows','originalLinearIndex', ...
    'selectedSamplePositions','pilotBaselineNRMSE9','pilotMeanGainClass', ...
    'pilotAnnualGainClass','pilotAnnualWinnerModel','pilotDesignWeight', ...
    'pilotStratumId','-v7.3');
fprintf('Pilot master saved: %s\n', cfg.pilotMasterFile);
fprintf('Assigned original grid rows: %d to %d (%d unique).\n', ...
    min(originalValidRows),max(originalValidRows),numel(unique(originalValidRows)));
T = table((1:cfg.nLocal)',originalValidRows,validLinearIndex,validLat,validLon, ...
    pilotStratumId,pilotMeanGainClass, 'VariableNames', ...
    {'PilotRow','OriginalValidRow','LinearIndex','Latitude','Longitude', ...
    'Stratum','MeanGainClass'});
writetable(T,fullfile(cfg.outputRoot,'pilot10_selected_cells.csv'));
disp(T(1:min(10,height(T)),:));
end
