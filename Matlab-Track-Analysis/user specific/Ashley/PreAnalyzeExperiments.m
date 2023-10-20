function [eset] = PreAnalyzeExperiments(eset)
%function [eset] = PreAnalyzeExperiments(eset)
%
%run some routines to check data and make sure good

%delete all tracks that start from the edge of the dish
%where second array is the limits that define "the edge", [xstart xend
%ystart yend]
eset=DeleteTracksFromEdge(eset,[400 1950 150 1750]);

%detect collisions
eset.executeExperimentFunction('detectCollisions', 30);

%delete all tracks that are too slow
eset=DeleteTracksThatAreSlow(eset);

%check for segmentation errors
[inds]=eset.executeExperimentFunction('detectPossibleSegmentationProblems');
for n=1:length([eset.expt])
    if (~isempty(inds{n}))
    disp(['Detected possible segmentation errors in expt ' num2str(n) ' and track ' num2str(inds{n}) '.']);
    end
end




end

