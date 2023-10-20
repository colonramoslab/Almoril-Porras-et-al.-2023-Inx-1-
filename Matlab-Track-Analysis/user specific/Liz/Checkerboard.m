%initialize variables
ecl = ESetCleaner;

%% Load experiment, if one already exists then don't load
    if (~exist('eset','var')) || (exist('esetReload','var') && esetReload)  % only load if you haven't already
        eset= ExperimentSet.fromFiles();
        esetReload=false; %if want to load another experiment, define esetReload = true;
    end
    
%% Load image of checkerboard used
if (~exist('checkerpics','var')) || (exist('checkerReload','var') && checkerReload)
    checkerpics = cell(1,length(eset.expt));
    for i = 1:length(eset.expt)
        [fn basedir] = uigetfile('*.tiff', 'Select the picture of the checkerboard used in each expt');
        filepath = strcat(basedir,fn);
        checkerpics{1,i} = imread(filepath);
    end
    checkerReload = false; %if want to load another array of pictures, define checkerReload = true;
end

% %% Create binary images of each checkerboard picture
% binarypics = cell(1,length(eset.expt));
% for i = 1:length(eset.expt)
%     binarypics{1,i} = checkerpics{1,i}>80;
% end
% 
% %% Define borders in each picture and store locations in array "borderlocs2"
% borderlocs2 = cell(1,length(eset.expt));
% for i = 1:length(eset.expt)
%     [r c] = size(checkerpics{1,i});
%     borderlocs2{1,i}= zeros(r,c);
%     for j= 2:r-1
%         for k = 2:c-1
%             if (binarypics{1,i}(j,k-1)~=binarypics{1,i}(j,k))||(binarypics{1,i}(j,k)~=binarypics{1,i}(j,k+1))
%                 borderlocs2{1,i}(j,k)=1; %%stores 1 in cells which are a vertical border
%             end
%             if (binarypics{1,i}(j-1,k)~=binarypics{1,i}(j,k))||(binarypics{1,i}(j,k)~=binarypics{1,i}(j+1,k))
%                 borderlocs2{1,i}(j,k)=2; %%stores 2 in cells which are a horizontal border
%             end
%         end
%     end
% end

%% Iterate through checkerpics and notice when borders are changing
%%different approach than borderlocs because it picks up gradual border
%%changes rather than relying on a single threshold
%%also different because it iterates checkerpics rather than binarypics

borderlocs = cell(1,length(eset.expt));
p = 20; %%number of pixels to evaluate across
thresh = 35; %%threshold change to accept
for i = 1:length(eset.expt)
    [r c] = size(checkerpics{1,i});
    borderlocs{1,i}= zeros(r,c);
    for j = p+1:r-p
        for k = p+1:c-p
            left = int16(checkerpics{1,i}(j,k-p));
            center = int16(checkerpics{1,i}(j,k));
            right = int16(checkerpics{1,i}(j,k+p));
            up = int16(checkerpics{1,i}(j-p,k));
            down = int16(checkerpics{1,i}(j+p,k));
            if (abs(left-center)>=thresh)||(abs(center-right)>=thresh)
                borderlocs{1,i}(j,k)=1; %%stores 1 in cells which are a vertical border
            end
            if (abs(up-center)>=thresh)||(abs(center-down)>=thresh)
                borderlocs{1,i}(j,k)=1; %%stores 1 in cells which are a horizontal border
            end
        end
    end
end

%% Iterate through borderlocs and mark the corners

side = 40; %%equal to the side length of the square at a corner
for i = 1:length(eset.expt)
    [r c] = size(checkerpics{1,i});
    for j = side+1:r-side
        for k = side+1:c-side
            left = borderlocs{1,i}(j,k-side);
            center = borderlocs{1,i}(j,k);
            right = borderlocs{1,i}(j,k+side);
            up = borderlocs{1,i}(j-side,k);
            down = borderlocs{1,i}(j+side,k);
            if ((left==1)&&(right==1)&&(up==1)&&(down==1))
                borderlocs{1,i}(j,k)=0; %%stores 2 in cells which are in the corner location
            end
        end
    end
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