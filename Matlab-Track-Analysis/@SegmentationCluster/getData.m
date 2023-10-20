function [data, mu, sigma2] = getData(sc, point_source, varargin)
%
%gets data based on data field and operations
%function data = getData(sc, point_source, varargin)
%
%outputs
%   DATA: a O x T matrix, where O is the number of observation types
%   (datafields) and T is the number of observations (points)
%   MU: an O x 1 vector of mean values for the DATA 
%   SIGMA2: an O x 1 matrix of standard deviations for the values of DATA
%inputs:
%   SC < SegmentationCluster
%   point_source: an array of tracks, or an experiment or an experimentset
%
%   optional arguments:
%   'stand', true/false (F) - if true, we subtract mu from data and
%       divide by sigma2 -- data(i,t) --> (data(i,t) - mu(i))/sigma2(i)
%               default: false
%   'mu', 'sigma2' : use these values for standardization
%   'cattracks', true/false (T) - if true, we return data as single matrix
%                             if false(default), each track is a cell
%   'knownPointsOnly', true/false (F) - only gather data from known points
%                           
%   'knownPoints' - use these known points instead of sc's
%
%   if point_source is empty, then knownPointsOnly is true, and we use
%   knownPoints.track as our point sources
sc = sc(1);
stand = false;
mu = [];
sigma2 = [];
cattracks = true;
knownPointsOnly = false;
knownPoints = sc.knownPoints;

varargin = assignApplicable(varargin);
existsAndDefault('point_source', []);

if (isempty(point_source))
    if (isempty(knownPoints))
        data = [];
        mu = [];
        sigma2 = [];
        return;
    end
    point_source = unique([knownPoints.track]);
    knownPointsOnly = true;
end

if (isa (point_source, 'ExperimentSet'))
    e = [point_source.expt];
    tr = [e.track];
end
if (isa (point_source, 'Experiment'))
    tr = [point_source.track];
end
if (isa (point_source, 'Track'))
    tr = point_source;
end

if (~exist('tr', 'var'))
    error ('point_source must be ExperimentSet(s), Experiment(s), or Track(s)');
end

if (knownPointsOnly)
    tr = intersect([knownPoints.track], tr);
end

data = cell(size(tr));
for j = 1:length(tr)
    if (knownPointsOnly)
        inds = [knownPoints([knownPoints.track] == tr(j)).inds];
        if (isempty(inds))
            continue;
        end
    else
        inds = 1:length(tr(j).getDerivedQuantity('eti'));
    end
    data{j} = zeros(length(sc.datafields), length(inds));
    for k = 1:length(sc.datafields)
        op = sc.operation{k};
        data{j}(k,:) = op(tr(j).getDerivedQuantity(sc.datafields{k}, false, inds));
    end
end
data2 = cell2mat(data);
if (stand || isempty(mu) || isempty(sigma2))
    [data2,mua,sigma2a] = standardize(data2);
    if (isempty(mu))
        mu = mua;
    end
    if (isempty(sigma2))
        sigma2 = sigma2a;
    end
end
if (cattracks)
    if (stand)
        data = data2;
    else
        data = cell2mat(data);
    end
else
    if (stand)
        for j = 1:length(data)
            data{j} = standardize(data{j}, mu, sigma2);
        end
    end
end