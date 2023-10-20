cc_location = 'e:\larvalco2\calibration images\cc_olfaction.mat';
if (~exist('cc_olfaction', 'var'))
    load(cc_location);
end
for j = 1:length(basedirs)
    eset(j) = ExperimentSet.loadTrimStitchAndSave(basedirs{j}, esetnames{j}, [], cc_olfaction);
end