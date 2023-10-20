setEthylAcetateFilesTemporal
inds = 5:7;
basedirs = basedirs(inds);
esetnames = esetnames(inds);
labelnames = labelnames(inds);
for j = 1:length(basedirs)
ea_eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
end
for j = 1:length(ea_eset)
ea_eset(j).executeTrackFunction('setSegmentSpeeds');
ea_eset(j).executeTrackFunction('segmentTrack');
end
tno.timerange = [60 1800];
tno.rampType = 'exponential';
tno.period = 600;
tno.fieldname = 'vocppm';
tno.timeBinSize = 30;
ts1 = tic;
tc = temporalCalculations(ea_eset, tno);
toc(ts1)