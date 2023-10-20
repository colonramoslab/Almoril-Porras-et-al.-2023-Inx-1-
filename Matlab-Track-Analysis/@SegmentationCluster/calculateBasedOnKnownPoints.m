function calculateBasedOnKnownPoints(sc)
%updates mean and covariance based on points known to be in the cluster
%function calculateBasedOnKnownPoints(sc)
%
%outputs: none (updates sc)
%inputs:
%   SC < SegmentationCluster

[cm,ccv] = sc.meanAndCovOfKnownPoints();
if (length(sc) > 1)
    for j = 1:length(sc)       
        sc(j).clustMean = cm{j};
        sc(j).clustCov = ccv{j};
    end
else
    sc.clustMean = cm;
    sc.clustCov = ccv;
end


