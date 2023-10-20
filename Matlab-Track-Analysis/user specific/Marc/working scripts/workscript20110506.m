setEthylAcetateFilesTemporal
if (~exist('eset', 'var'))
    eset = ExperimentSet.fromMatFiles(fullfile(basedirs{3},'matfiles', esetnames{3}));
    eset.executeTrackFunction('segmentTrack');
    eset.expt.addTonToff('vocppm', 'triangle');
end

