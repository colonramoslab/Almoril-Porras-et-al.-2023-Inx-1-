function addKnownPoints(sc, track, inds)
%adds points to the list of those known to be part of the cluster
%function addKnownPoints(sc, track, inds)
%
%ouputs: none (modifies sc)
%inputs:
%   sc < SegmentationCluster
%   track < Track
%   inds -- indices (interpolated) of known points

kp.track = track;
kp.inds = inds;
sc.knownPoints = [sc.knownPoints kp];
