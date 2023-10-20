eset = ExperimentSet.fromFiles;

maxFrameInterval=4;
maxDist=8;
eset.executeExperimentFunction('stitchTracks', maxFrameInterval, maxDist, 'interactive', false);

ecl = ESetCleaner();
ecl.minSpeed = 0.4; %minimum average speed of the track to be kept
ecl.minDist = 80; %minimum distance must travel from the starting point to be kept
ecl.minPts = 400; %minimum number of points to be kept
ecl.clean(eset);

eset.executeExperimentFunction('segmentTracks');

starttimerange = [0 1800];
startrect = [700 300 1900 1600];
eset.executeExperimentFunction ('pruneTracks', starttimerange, startrect);

figure
hold on;
eset.executeTrackFunction('plotPath', 'displacement'); % plot all paths
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
tr=eset.expt(1).track;
tr(1).plotPath('displacement','r');
tr=eset.expt(6).track;
tr(1).plotPath('displacement','g');
axis equal;


figure;
hold on;
num=20;
for i=1:length(eset.expt);
    for j=1:length(eset.expt(i).track)
        if eset.expt(i).track(j).npts>3610 && eset.expt(i).track(j).startFrame<200 && num>0
            tr=eset.expt(i).track(j);
            if num==3
                tr.plotPath('displacement','r','inds',1:3600);
            else
                tr.plotPath('displacement','b','inds',1:3600);
            end
            num=num-1;
            hold on;
        end
    end
end
axis equal;
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0

        
