%Script to analyze larval tracks navigating a spatial light gradient
%produced by a digital projector. 

% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and temperature (.tmp) files, is loaded into
% MATLAB, then processed to show how larval behavior depends on 
% temperature.

%initialize variables
ecl = ESetCleaner;

% %load more than one eset
% esetNum=input('Enter number of esets you wish to load:')

%% Load experiment, if one already exists then don't load
    if (~exist('eset','var')) || (exist('esetReload','var') && esetReload)  % only load if you haven't already
        eset= ExperimentSet.fromFiles();
        esetReload=false; %if want to load another experiment, define esetReload = true;
    end
%% Clean up tracks and segment tracks
% %the below code can be used to display the eset cleaning levels 
% %graphically and to change default values specified below
% ecl.getReport(eset);

%throw out first minute of data acquistion
eset.executeExperimentFunction('trimTracks',[60 1200],[]);

%initializes variables for eset cleaning and cleans eset
if (~exist('cleanEset','var')) || (exist('reclean','var') && reclean)   
    ecl.minHTValid= .95;  %allow for user value inputs?
    ecl.minDist=75;
    ecl.minSpeed = 0.75;
    ecl.minPts = 500;
    ecl.clean(eset);
    reclean=false;
    cleanEset=false;
end

%fixes head/tail orientation
if (~exist('fixht','var')|| (exist('fixht','var') && fixht))
    eset.executeTrackFunction('fixHTOrientation');
    fixht=false;    
end

% Set segmentation threshold speeds )
if (~exist('autosetspeeds','var') || (exist('autosetspeeds','var') && autosetspeeds))
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% Segment the tracks
if (~exist('segmentTest','var') || (exist('segmentTest','var') && segmentTest))
    eset.executeTrackFunction('segmentTrack');
    segmentTest = false; 
end

%% Plot some pretty graphs

%initialize fignum to zero
fignum = 0;

%plot all paths for eset
fignum = fignum + 1; figure(fignum); clf(fignum);
eset.executeTrackFunction('plotPath');

% 
% 
% fignum = fignum + 1; figure(fignum); clf(fignum);
% eset.makeReorientationHistogram('theta',deg2rad(-180:30:180),'r2d', true);
% 
% fignum = fignum + 1; figure(fignum); clf(fignum);
% eset.makeReorientationHistogram('theta',deg2rad(-180:30:180),'minHS', 1, 'r2d', true);
% title
% 
% hd=eset.gatherSubField('headSwing','headDir');
% td=eset.gatherSubField('headSwing','tailDir');
% dt = eset.gatherSubField('headSwing', 'maxTheta');
% 
% acc=logical(eset.gatherSubField('headSwing','accepted'));
% 
% startUp=(abs(td-pi/2)<pi/4);
% startDown=(abs(td+pi/2)<pi/4);
% startLeft=(abs(td)<pi/4);
% startRight=(abs(td)>3*pi/4);
% 
% turnLeft = dt > 0;
% turnRight = dt < 0;
% 
% upRightRate = mean(acc(startUp & turnRight));
% upLeftRate = mean(acc(startUp & turnLeft));
% downRightRate= mean(acc(startDown & turnRight));
% downLeftRate= mean(acc(startDown & turnLeft));
% leftRightRate = mean(acc(startLeft & turnRight));
% leftLeftRate = mean(acc(startLeft & turnLeft));
% rightLeftRate = mean(acc(startRight & turnLeft));
% rightRightRate = mean(acc(startRight & turnRight));
% 
% %make bar graph
% fignum = fignum + 1; figure(fignum); clf(fignum);
% 
% bar([[upRightRate;upLeftRate],[downRightRate;downLeftRate],[leftRightRate;leftLeftRate],[rightRightRate;rightLeftRate]]');
% set(gca,'XTickLabel',{'up','down','left','right'})
% legend('right','left');
% title ('probability of accepting a head sweep');
% ylabel ('probability of acceptance');
% 
% fignum = fignum + 1; figure(fignum); clf(fignum);
% 
% bar([[sum(startUp&turnRight);sum(startUp&turnLeft)],[sum(startDown&turnRight);sum(startDown&turnLeft)]]');
% set(gca,'XTickLabel',{'up','down','left','right'})
% legend('right','left');
% title ('total number of head sweeps');
% 
% td=eset.gatherSubField('firstHS','tailDir');
% dt = eset.gatherSubField('firstHS', 'maxTheta');
% 
% acc=logical(eset.gatherSubField('firstHS','accepted'));
% 
% 
% startUp=(abs(td-pi/2)<pi/4);
% startDown=(abs(td+pi/2)<pi/4);
% startLeft=(abs(td)<pi/4);
% startRight=(abs(td)>3*pi/4);
% 
% turnLeft = dt > 0;
% turnRight = dt < 0;
% 
% fignum = fignum + 1; figure(fignum); clf(fignum);
% 
% bar([[sum(startUp&turnRight);sum(startUp&turnLeft)],[sum(startDown&turnRight);sum(startDown&turnLeft)],...
%     [sum(startLeft&turnRight);sum(startLeft&turnLeft)],[sum(startRight&turnRight);sum(startRight&turnLeft)]]');
% set(gca,'XTickLabel',{'up','down','left','right'})
% legend('right','left');
% title ('total number of first head sweeps');
% %%
% %plot(td(acc), hd(acc), 'g.',td(~acc),hd(~acc),'r.');
