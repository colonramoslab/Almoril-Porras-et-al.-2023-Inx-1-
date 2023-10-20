function cases = createBnetCases (sm, tracks, varargin)
%creates cases for training bayesnet hidden markov model
%function cases = createBnetCases(sm, tracks, varargin)
%
%outputs:
%   CASES - a 1 x length(tracks) cell array, cases{i} is
%           a cell array of size (3 x num interped pts in track)
%           cases(1,j) = seg cluster # if known (marked by user) at
%               interped pt j
%           cases(2,j) = [] - in a mixture model this is the gaussian
%               picked from -- never known
%           cases(3,j) = data (always observed) from interped pt j
%           
%inputs:
%   SM < SegmentationModel
%   tracks < Track
%
%optional arguments:
%   none

cases = cell(size(tracks));

for j = 1:length(tracks)
    cases{j} = makeCase(sm, tracks(j), varargin{:});
end

function cs = makeCase(sm, track, varargin)
cs = cell(3, length(track.getDerivedQuantity('eti')));
sc = sm.segmentationClusters;
data = zeros(length(sc(1).datafields), length(track.getDerivedQuantity('eti')));
for j = 1:length(sc(1).datafields)
    op = sc(1).operation{j};
    data(j,:) = op(track.getDerivedQuantity(sc(1).datafields{j}));
end

for j = 1:length(sc)
    kp = sc(j).knownPoints([sc(j).knownPoints.track] == track);
    if (~isempty(kp))
        size([kp.inds])
        [cs{1, [kp.inds]}] = deal(j);
    end
end
for j = 1:size(data,2)
    cs{3,j} = data(:,j);
end