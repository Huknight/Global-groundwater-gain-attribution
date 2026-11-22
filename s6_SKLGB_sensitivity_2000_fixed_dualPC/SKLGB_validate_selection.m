function SKLGB_validate_selection(sampleFile)
% Lightweight preflight. Checks nested identities and 18-stratum counts.
P=load(sampleFile,'selectedValidRow','selectedLinearIndex', ...
    'meanGainClass','stratumId','populationStratumCounts');
N=numel(P.selectedValidRow);
rows=double(P.selectedValidRow(:));
ix=double(P.selectedLinearIndex(:));
if numel(unique(rows))~=N || numel(unique(ix))~=N || ...
        any(diff(rows)<=0) || any(diff(ix)<=0)
    error('Sample file has duplicate or unordered original-grid indices.');
end
if N~=400 && N~=2000
    error('Expected a 400- or 2000-cell sample MAT; got %d.',N);
end
if any(P.meanGainClass(:)<1 | P.meanGainClass(:)>6) || ...
        any(P.stratumId(:)<1 | P.stratumId(:)>18)
    error('Invalid gain class or stratum identifiers.');
end
if any(double(P.stratumId(:)) ~= ...
        3*(double(P.meanGainClass(:))-1) + ...
        mod(double(P.stratumId(:))-1,3) + 1)
    error('Invalid class-stratum encoding.');
end

if N==2000
    priorityInfo=load(sampleFile,'nestedPriority');
    priority=double(priorityInfo.nestedPriority(:));
    if ~isequal(sort(priority),(1:N)')
        error('nestedPriority must be a permutation of 1:2000.');
    end
    tenRows=rows(priority<=10);
    fourRows=rows(priority<=400);
    expectedOriginalTen=[5;2022;3806;5692;7613;9576;11175;13233;15510;17275];
    if ~isequal(tenRows(:),expectedOriginalTen)
        error('The 10-grid pilot was not preserved.');
    end
    old400File=fullfile(fileparts(sampleFile), ...
        'SKLGB_hyperparameter_sample400.mat');
    if isfile(old400File)
        A=load(old400File,'selectedValidRow');
        if ~isequal(fourRows(:),double(A.selectedValidRow(:)))
            error('The original 400-grid sample was not preserved.');
        end
    end
    fprintf('Nested sampling validated: original 10 contained in 400 contained in 2000.\n');
end
fprintf('Sampling preflight: %d distinct original grid cells; %d strata observed.\n', ...
    N,numel(unique(P.stratumId(:))));
end
