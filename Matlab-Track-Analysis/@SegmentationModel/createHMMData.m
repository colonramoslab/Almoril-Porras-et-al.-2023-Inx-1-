function [data, known_points] = createHMMData (sm, tracks, varargin)
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
%   'stand', true/false (T) - if true, we subtract mu from data and
%       divide by sigma2 -- data(i,t) --> (data(i,t) - mu(i))/sigma2(i)
%               default: true
%    if standardize is true, we use hmm_musub and hmm_signorm to
%    standardize
%    if these are unavailable, we do not standardize

stand = true;
varargin = assignApplicable(varargin);
sc = sm.segmentationClusters;
known_points = cell(size(tracks));
if (isempty(sm.hmm_musub) || isempty(sm.hmm_sigmnorm) || stand == false)
    data = sc(1).getData(tracks, 'cattracks', false, 'stand', false);
else
    data = sc(1).getData(tracks, 'cattracks', false, 'stand', true, 'mu', sm.hmm_musub, 'sigma2', sm.hmm_sigmnorm);
end

for k = 1:length(tracks)
    known_points{k} = zeros(1, length(tracks(k).getDerivedQuantity('eti')));
    for j = 1:length(sc)
        kp = sc(j).knownPoints([sc(j).knownPoints.track] == tracks(k));
        
        if (~isempty(kp))

            known_points{k}([kp.inds]) = j;
        end
    end
end


%{
data = cell(size(tracks));
known_points = cell(size(tracks));
for j = 1:length(tracks)
    [d, kp] = makeHMMData(sm, tracks(j));
    data{j} = d;
    known_points{j} = kp;
end

function [data, known_points] = makeHMMData(sm, track, varargin)

sc = sm.segmentationClusters;
data = zeros(length(sc(1).datafields), length(track.getDerivedQuantity('eti')));
for j = 1:length(sc(1).datafields)
    op = sc(1).operation{j};
    data(j,:) = op(track.getDerivedQuantity(sc(1).datafields{j}));
end
known_points = false(1, length(track.getDerivedQuantity('eti')));
for j = 1:length(sc)
    kp = sc(j).knownPoints([sc(j).knownPoints.track] == track);
    if (~isempty(kp))
        known_points(kp.inds) = true;
    end
end
%}