%first example, loading a single experiment
mintracklength = 500;

%load experiment set, dont load if one already has been loaded
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles('\\labnas1\Share\Ashley Extracted\Checkerboard\main experiment','minpts',mintracklength);
    %use the below code to get report when analyzing data for the first
    %time
    ecl = ESetCleaner;
    ecl.getReport(eset);
end


%clean up bad tracks and treshold data
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid=.95;
    ecl.minDist=0;
    ecl.minSpeed = 2;
    ecl.clean(eset);
    
    trimRect = [355 800 2145 1795];
    eset.executeExperimentFunction('trimTracks', [], trimRect);
    
    cleanEset = false;
    
end

%ecl.getReport(eset); can be used to get report to obtain thresh values
existsAndDefault('fixht','true');
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');  
    fixht = false;
end

existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end
%if there isn't a variable called segment, create a variable called segment
%and set it to true
existsAndDefault('segment', true);

%set segmentation speeds automatically and segment tracks
if (segment)
    eset.executeTrackFunction('segmentTrack');
    segment = false; %don't do it again when you run the script next
end
    
%inds = eset.executeExperimentFunction('detectPossibleSegmentationProblems')

fignum = 0;

fignum = fignum + 1; figure(fignum); clf(fignum);
eset.makeReorientationHistogram('theta',deg2rad(-180:30:180),'r2d', true);

fignum = fignum + 1; figure(fignum); clf(fignum);
eset.makeReorientationHistogram('theta',deg2rad(-180:30:180),'minHS', 1, 'r2d', true);

hd=eset.gatherSubField('headSwing','headDir');
td=eset.gatherSubField('headSwing','tailDir');
dt = eset.gatherSubField('headSwing', 'maxTheta');

acc=logical(eset.gatherSubField('headSwing','accepted'));

startUp=(abs(td-pi/2)<pi/4);
startDown=(abs(td+pi/2)<pi/4);
startLeft=(abs(td)<pi/4);
startRight=(abs(td)>3*pi/4);

turnLeft = dt > 0;
turnRight = dt < 0;

upRightRate = mean(acc(startUp & turnRight));
upLeftRate = mean(acc(startUp & turnLeft));
downRightRate= mean(acc(startDown & turnRight));
downLeftRate= mean(acc(startDown & turnLeft));
leftRightRate = mean(acc(startLeft & turnRight));
leftLeftRate = mean(acc(startLeft & turnLeft));
rightLeftRate = mean(acc(startRight & turnLeft));
rightRightRate = mean(acc(startRight & turnRight));

%make bar graph
fignum = fignum + 1; figure(fignum); clf(fignum);

bar([[upRightRate;upLeftRate],[downRightRate;downLeftRate],[leftRightRate;leftLeftRate],[rightRightRate;rightLeftRate]]');
set(gca,'XTickLabel',{'up','down','left','right'})
legend('right','left');
title ('probability of accepting a head sweep');
ylabel ('probability of acceptance');

fignum = fignum + 1; figure(fignum); clf(fignum);

bar([[sum(startUp&turnRight);sum(startUp&turnLeft)],[sum(startDown&turnRight);sum(startDown&turnLeft)]]');
set(gca,'XTickLabel',{'up','down','left','right'})
legend('right','left');
title ('total number of head sweeps');

td=eset.gatherSubField('firstHS','tailDir');
dt = eset.gatherSubField('firstHS', 'maxTheta');

acc=logical(eset.gatherSubField('firstHS','accepted'));


startUp=(abs(td-pi/2)<pi/4);
startDown=(abs(td+pi/2)<pi/4);
startLeft=(abs(td)<pi/4);
startRight=(abs(td)>3*pi/4);

turnLeft = dt > 0;
turnRight = dt < 0;

fignum = fignum + 1; figure(fignum); clf(fignum);

bar([[sum(startUp&turnRight);sum(startUp&turnLeft)],[sum(startDown&turnRight);sum(startDown&turnLeft)],...
    [sum(startLeft&turnRight);sum(startLeft&turnLeft)],[sum(startRight&turnRight);sum(startRight&turnLeft)]]');
set(gca,'XTickLabel',{'up','down','left','right'})
legend('right','left');
title ('total number of first head sweeps');
%%
%plot(td(acc), hd(acc), 'g.',td(~acc),hd(~acc),'r.');

