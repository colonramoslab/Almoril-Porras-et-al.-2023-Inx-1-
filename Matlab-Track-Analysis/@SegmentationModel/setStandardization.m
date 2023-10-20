function setStandardization(sm, tracklist)
if (isempty(tracklist))
    kp = [sm.segmentationClusters.knownPoints];
    tracklist = unique([kp.track]);
end
[~,sm.hmm_musub, sm.hmm_sigmnorm] = sm.segmentationClusters(1).getData(tracklist, 'stand', true);
