function [cm, ccv, cmr] = meanAndCovOfKnownPoints(sc, numclusters, varargin)
%calculates mean and covariance based on points known to be in the cluster
%function calculateBasedOnKnownPoints(sc, numclusters)
%
%outputs: 
%   CM: mean of data for known points
%   CCV: covariance of data for known points
%   cmr: mix ratio of clusters
%   if sc is an array, CM,CCV are cells with the same length as sc
%inputs:
%   SC < SegmentationCluster
%   numclusters (defaults to 1): if > 1, fits to a gaussian mixture model
%   optional arguments:
%       'stand' - if true standardize data
%       'mu', 'sigma2' - use these values to standardize (recommended to
%       pass)

existsAndDefault('numclusters', 1);

if (length(sc) > 1)
    for j = 1:length(sc)
        [out1,out2,out3] = meanAndCovOfKnownPoints(sc(j), numclusters, varargin{:});
        cm{j} = out1;
        ccv{j} = out2;
        cmr{j} = out3;
    end
    return
end

%{
kp = sc.knownPoints;
data = zeros(length(sc.datafields), length([kp.inds]));

if isempty(sc.operation)
    sc.operation = {@(x) x};
    for j = 1:length(sc.datafields)
        sc.operation{j} = @(x) x;
    end
end

ind0 = 0;
for k = 1:length(kp)
    inds = ind0 + (1:length(kp(k).inds));
    ind0 = max(inds);
    for j = 1:length(sc.datafields)
        op = sc.operation{j};
        data(j,inds) = op(kp(k).track.getDerivedQuantity(sc.datafields{j},false,kp(k).inds));
    end
end
%}

data = sc.getData([], varargin{:});
gmfit = gmdistribution.fit(data', numclusters,'Regularize',1E-3);
cmr = gmfit.PComponents;
ccv = gmfit.Sigma;
cm = gmfit.mu;
  
