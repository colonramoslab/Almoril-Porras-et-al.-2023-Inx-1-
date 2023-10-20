function weightedEM(sm, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ptsource = sm.eset;
sc = sm.segmentationClusters;
maxsteps = 100;
varagin = assignApplicable(varargin);
termeps = 1E-4;
%find known cluster means and covariances; these do not change 
sc.calculateBasedOnKnownPoints;
[kcm, kccv] = sc.meanAndCovOfKnownPoints();
oldpp = zeros(size(sc.posteriorProbabilities(ptsource)));
delta = zeros([1 maxsteps]);
for j = 1:maxsteps
    %match points to clusters, based on previous guess
    pp = sc.posteriorProbabilities(ptsource);
    priorprob = sum(pp,2)/size(pp,2);
    delta(j) = sum(abs(pp(:)-oldpp(:)))/numel(pp);
    if  delta(j) < termeps
        break;
    end
    oldpp = pp;
    %update guess based on matches
    [gcm, gccv] = sc.meanAndCovEMStep(ptsource, pp);    
    for k = 1:length(kcm)
        sc(k).priorProbability = priorprob(k);
        sc(k).clustMean = sm.weightOnKnown * kcm{k} + (1-sm.weightOnKnown)*gcm{k};
        sc(k).clustCov = sm.weightOnKnown * kccv{k} + (1-sm.weightOnKnown)*gccv{k};
    end
end
plot (delta)