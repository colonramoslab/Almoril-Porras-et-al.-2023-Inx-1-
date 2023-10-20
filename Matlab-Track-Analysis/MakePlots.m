function [time, x, y]=MakePlots

eset = ExperimentSet.fromFiles;

maxFrameInterval=3;
maxDist=4.5;
eset.executeExperimentFunction('stitchTracks', maxFrameInterval, maxDist, 'interactive', false);

ecl = ESetCleaner();
ecl.minSpeed = 0.4; %minimum average speed of the track to be kept
ecl.minDist = 80; %minimum distance must travel from the starting point to be kept
ecl.minPts = 400; %minimum number of points to be kept
ecl.clean(eset);


timerange = [0 1800];
validrect = [1000 25 2000 1850];
% timerange = [];
%validrect = [1000 500 2000 1500];
eset.executeExperimentFunction('trimTracks', timerange, validrect);

eset.executeExperimentFunction('segmentTracks');


timeAxis = 0:100:1800;
[time, centerofmass] = eset.meanField2vsField1('eti', 'sloc', timeAxis);
xposition = centerofmass(1,:);
x=xposition - xposition(1,1);
yposition = centerofmass(2,:);
y=yposition - yposition(1,1);

figure(1);
hold on;
eset.executeTrackFunction('plotPath');
hold off;
axis equal;

figure(2);
hold on;
eset.executeTrackFunction('plotPath', 'displacement'); % plot all paths
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
hold off;
axis equal;

figure(3);
rundirection = eset.gatherSubField('run', 'meanTheta');
runduration = eset.gatherSubField('run', 'runTime');
plot (rad2deg(rundirection), runduration,'.');

figure(4);
thetaAxis = deg2rad(-180:30:180);
eset.makeHistogram('theta',thetaAxis,'runs','polar',true,'r2d',true);

figure(5);
%to exclude first 5 and last 5 pts from run
runtheta = eset.gatherFromSubField('run', 'theta', 'trimpts', 5);
runspeed = eset.gatherFromSubField('run', 'speed', 'trimpts', 5);
[newtheta, meanspeed, stderrorspeed] = meanyvsx(runtheta, runspeed, thetaAxis);
errorbar(rad2deg(newtheta), meanspeed, stderrorspeed);

figure(6);
hold on
runstart = eset.gatherSubField('run', 'startTheta');
runend = eset.gatherSubField('run', 'endTheta');
rundt = diff(unwrap([runstart;runend])); % this takes care of 2 pi discontinuities
dthetaAxis = deg2rad(-180:15:180);
startedRight = abs(runstart) < pi / 4;  
startedLeft=runstart>3*pi/4 | runstart<-3*pi/4;
startedUp=runstart>pi/4 & runstart<3*pi/4;
startedDown=runstart>-3*pi/4 & runstart<-pi/4;
plot (rad2deg(dthetaAxis), hist(rundt(startedRight), dthetaAxis)./sum(hist(rundt(startedRight), dthetaAxis)),'r');
plot (rad2deg(dthetaAxis), hist(rundt(startedLeft), dthetaAxis)./sum(hist(rundt(startedLeft), dthetaAxis)),'b');
plot (rad2deg(dthetaAxis), hist(rundt(startedUp), dthetaAxis)./sum(hist(rundt(startedUp), dthetaAxis)),'m');
plot (rad2deg(dthetaAxis), hist(rundt(startedDown), dthetaAxis)./sum(hist(rundt(startedDown), dthetaAxis)),'g');
hold off

run1=[];
run2=[];
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track);
        for k=2:length(eset.expt(i).track(j).run);
            run1=[run1, eset.expt(i).track(j).run(k-1).endTheta];
            run2=[run2, eset.expt(i).track(j).run(k).startTheta];
        end
    end
end
dthetaAxis = deg2rad(-180:15:180);
rundt = diff(unwrap([run1;run2]));
startedRight = abs(run1) < pi / 4;  
startedLeft=run1>3*pi/4 | run1<-3*pi/4;
startedUp=run1>pi/4 & run1<3*pi/4;
startedDown=run1>-3*pi/4 & run1<-pi/4;

figure(7)
hold on;
plot (rad2deg(dthetaAxis), hist(rundt(startedRight), dthetaAxis)./sum(hist(rundt(startedRight), dthetaAxis)),'r');
plot (rad2deg(dthetaAxis), hist(rundt(startedLeft), dthetaAxis)./sum(hist(rundt(startedLeft), dthetaAxis)),'b');
plot (rad2deg(dthetaAxis), hist(rundt(startedUp), dthetaAxis)./sum(hist(rundt(startedUp), dthetaAxis)),'m');
plot (rad2deg(dthetaAxis), hist(rundt(startedDown), dthetaAxis)./sum(hist(rundt(startedDown), dthetaAxis)),'g');
hold off



