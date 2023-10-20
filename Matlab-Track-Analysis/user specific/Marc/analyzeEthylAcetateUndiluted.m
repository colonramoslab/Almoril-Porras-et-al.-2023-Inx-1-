basedir = 'D:\Marc Processed\maggots\ethyl acetate undil\';
d = dir([basedir '*.bin']);
maxToLoad = 6; %allows us to load just a subset for testing

if (~exist('expt','var'))
    for j = 1:(min(length(d),maxToLoad))
        n = [basedir d(j).name];
        s = strfind(n,'.');
        n2 = [n(1:s(end)) 'tim'];
        expt(j) = Experiment.fromFile(n, n2, true, [], 100);
    end
    
end
if (~exist('eset', 'var'))
    eset = ExperimentSet();
    eset.expt = expt;
    maxFrameInterval = 5;
    maxDist = 15;
    eset.executeExperimentFunction ('stitchTracks', maxFrameInterval, maxDist);
    eset.executeTrackFunction('fixHTOrientation');
    %analysis rect = 300 100 2275 1940
    validrect = [325 125 2250 1915];
    %discard any tracks that start at the edges
    eset.executeExperimentFunction ('pruneTracks', [], validrect); 
    %discard first 60 seconds of the experiment
    eset.executeExperimentFunction ('trimTracks', [60 Inf], []);
end
t = [eset.expt.track];
meanspeed = zeros(size(t));
for j = 1:length(t)
    meanspeed(j) = mean(t(j).getDerivedQuantity('speed'));
end
%figure(1); clf(1)
%hist (meanspeed,0:0.1:7);

figure(1); clf(1);
sx = 0:0.1:10;
npts = length(eset.gatherField('speed'));
bar (sx, eset.makeHistogram('speed', sx)/npts);
so = expt(1).so;
hold all

%plot (sx, hist(meanspeed, sx)/length(meanspeed), 'k-', 'LineWidth', 3);

for j = 1:length(eset.expt)
    plot (sx, hist(eset.expt(j).gatherField('speed'), sx)/length(eset.expt(j).gatherField('speed')),'LineWidth',3);
end
legend('all', '1', '2', '3', '4', '5');

yl = get(gca, 'YLim');
plot (so.stop_speed_cut * [1 1], yl, 'r--', so.start_speed_cut * [1 1], yl, 'g--'); 
hold off


minpts = 400;
for j = 1:length(eset.expt)
    eset.expt(j).track = eset.expt(j).track([eset.expt(j).track.npts] > minpts);
end



if (isempty([eset.expt(1).track.run]) || exist('resegment','var') && resegment)    
    for j = 1:length(eset.expt)
        so = eset.expt(j).so;
        pctl = percentile(eset.expt(j).gatherField('speed'), [0.15 0.25]);
        so.start_speed_cut = pctl(2);
        so.stop_speed_cut = pctl(1);
        eset.expt(j).so = so;
        eset.executeTrackFunction('segmentTrack', so);
    end
    resegment = false;
end

tx = deg2rad(-180:45:180);
thetadist = eset.makeHistogram('theta', tx, 'runs');
thetadist(1) = thetadist(1)+thetadist(end);
thetadist(end) = thetadist(1);
%thetadist = thetadist(1:end-1);
thetaend = eset.makeSubFieldHistogram('run','endTheta', tx);
thetaend(1) = thetaend(1) + thetaend(end);
thetaend(end) = thetaend(1);
%thetaend = thetaend(1:end-1);
%tx = tx(1:end-1)
figure(2); clf(2)
bar (rad2deg(tx), thetaend./thetadist * 240); title 'Instantaneous Reorientation Rate vs. Heading';
ylabel('rate (min^{-1})'); xlabel ('heading (-90 is towards higher [odor])'); embiggen();

figure(3); clf(3);
sthist = eset.makeSubFieldHistogram('run', 'startTheta', tx); 
sthist(1) = sthist(1)+sthist(end); sthist(end) = sthist(1);
bar (rad2deg(tx), sthist); title ('Initial Run Direction Histogram');
ylabel ('num runs'); xlabel ('heading (-90 is towards higher [odor])'); embiggen();
title ('histogram of initial run headings')

figure(4);
sthist = eset.makeSubFieldHistogram('run', 'meanTheta', tx); 
sthist(1) = sthist(1)+sthist(end); sthist(end) = sthist(1);
bar (rad2deg(tx), sthist); title ('Mean Run Direction Histogram');
ylabel ('num runs'); xlabel ('heading (-90 is towards higher [odor])'); embiggen();
