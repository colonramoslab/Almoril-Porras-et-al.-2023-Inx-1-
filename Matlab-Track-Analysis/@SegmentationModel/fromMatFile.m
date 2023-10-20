function sm = fromMatFile(fn, eset)
%function sm = fromMatFile(fn, eset)

existsAndDefault('eset', []);

load (fn, '-mat', 'segmentationModel');

if (~existsAndDefault('segmentationModel', segmentationModel))
    warn (['no variable called segmentationModel found in ' fn]);
    sm = segmentationModel;
    return;
end

sm = segmentationModel.fromSaveFriendly(eset);
