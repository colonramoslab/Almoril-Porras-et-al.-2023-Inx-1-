%Script to analyze larval tracks navigating a temporal light gradient
%produced by a digital projector. 

% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and light files, is loaded into
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
%eset.executeExperimentFunction('trimTracks',[60 1200],[]);

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