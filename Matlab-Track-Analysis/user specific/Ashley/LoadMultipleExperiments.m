function [eset] = LoadMultipleExperiments ()
%function [eset] = LoadMultipleExperiments ()
%
%Load a bunch of tracks.bin into the experiment. Also fix the head and
%tails.

eset=ExperimentSet();

[FileName, DataPath] = uigetfile('*.*','Select the .bin files you would like to use', 'MultiSelect', 'on');

numfiles=iscellstr(FileName)*length(FileName);
if (numfiles==0) numfiles=1; end

for n=1:numfiles
    disp(['Loading file #' num2str(n)]);
    if (numfiles==1)
        fn=[DataPath FileName];
    else
        fn=[DataPath FileName{n}];
    end
    timfn=[fn(1:end-4) '.tim'];
    expt(n) = Experiment.fromFile (fn, timfn, true, [], 400);
    expt(n).executeTrackFunction('fixHTOrientation');
    expt(n).executeTrackFunction('setSegmentSpeeds');
    expt(n).executeTrackFunction('segmentTrack');
    expt(n).executeTrackFunction('calculateDerivedQuantity', {'smid'})
    expt(n).executeTrackFunction('calculateDerivedQuantity', {'vmid'})
end

eset.expt=expt;


end