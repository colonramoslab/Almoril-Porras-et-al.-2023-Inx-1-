function toMatFile (sm, filename)
%function toMatFile (sm, filename)

segmentationModel = sm.toSaveFriendly();
save (filename, 'segmentationModel', '-v6'); %no compression, takes forever to load
