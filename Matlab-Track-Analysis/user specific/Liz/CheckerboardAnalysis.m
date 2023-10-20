%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CHECKERBOARD ANALYSIS SCRIPT
%
%   Written by Liz    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to analyze larval tracks navigating a checkerboard light stimulus
% made by a digital projector.

% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and light (.bin) files, is loaded into
% MATLAB

%created resetAll as a var to reload a new experiment from scratch

%% LOAD EXPERIMENT
    if (~exist('eset','var')) || (exist('esetReload','var') && esetReload) || exist('resetAll','var')  % only load if you haven't already
        eset= ExperimentSet.fromFiles();
        esetReload=false; %if want to load another experiment, define esetReload = true;
    end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP AND SAVE TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %the below code can be used to display the eset cleaning levels 
% %graphically and to change default values specified below
% ecl.getReport(eset);

% %Only look at first 23 minutes of data
% %have to do this before you do anything 
% if (~exist('trim','var') || (exist('trim','var') && trim)) || exist('resetAll','var')
%     eset.executeExperimentFunction('trimTracks',[0 endSec],[]);
%     trim = false; 
% end

%init variables
ecl = ESetCleaner;

if (~exist('reclean','var')) || (exist('reclean','var') && reclean) || exist('resetAll','var') %clean if you haven't already
    ecl.minHTValid= .95; 
    ecl.minDist=75;
    ecl.minSpeed = 0.75;
    ecl.minPts = 500;
    ecl.clean(eset);
    reclean=false;  %reset to true if you want to reclean
end

%close figures opened by ecl.clean(eset)
close all

%fixes head/tail orientation
if (~exist('fixht','var')|| (exist('fixht','var') && fixht)) || exist('resetAll','var') %fix ht if you haven't already
    eset.executeTrackFunction('fixHTOrientation');
    fixht=false;    %reset to true if you want to refix h/t orientation
end

%set headsweep theta min to 15 deg for segmentation (found this doesn't
%miss as many headsweeps for second instars
for j=1:length(eset)
    eset.expt(j).so.headswing_start=deg2rad(15);
end

% Set segmentation threshold speeds
if (~exist('setSpeeds','var') || (exist('setSpeeds','var') && setSpeeds)) || exist('resetAll','var')
    eset.executeTrackFunction('setSegmentSpeeds');
    setSpeeds = false;
end

% Save to mat file

% % Segment the tracks
% if (~exist('segmentTest','var') || (exist('segmentTest','var') && segmentTest)) || exist('resetAll','var')
%     eset.executeTrackFunction('segmentTrack');
%     segmentTest = false; 
% end
%     
% Need to check segmentation by playing some movies here
