%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRECTIONAL LIGHT ANALYSIS SCRIPT
%
%   Written by Liz    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to analyze larval tracks navigating a directional light gradient
% produced by a projector

% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and light (.bin) files, is loaded into
% MATLAB

endSec = 1400; %cut data 23 minutes (derived by looking at all time and deciding that response was max during first 23 mins)
trimBuffer=50; %pixel buffer distance to trim within extraction window

%create resetAll as a var to reload a new experiment from scratch

%% LOAD EXPERIMENT
    if (~exist('eset','var')) || exist('resetAll','var')  % only load if you haven't already
        matF=input('Do you want to load from MAT Files. Enter 1 or 0: ');
        if matF
            eset=LoadFromMat();
            bypass=true;
        else
            eset = ExperimentSet.fromFiles();
            bypass=false;
        end
    end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP AND SEGMENT TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %the below code can be used to display the eset cleaning levels 
% %graphically and to change default values specified below
% ecl.getReport(eset);

% %throw out first minute of data acquistion
% eset.executeExperimentFunction('trimTracks',[60 1200],[]);

%Only look at first 1400 secs of data - determined this was the range
%where larvae remain maximally sensitive to light changes
%have to do this before you do anything 
%trim a little past 1400 so that reorientation plots for last time bin
%aren't erroneous

%analysis options struct keeps track of which analyses have been run on the
%data you have open

if ~exist('ao','var')
    ao=[];
end

if ~bypass
   
    if ~isfield(ao,'trim') || ao.trim==false
        disp(['trimming track time to be less than ' num2str(endSec) ' and smaller than extraction window']);
        loc=eset.gatherSubField('pt','loc'); %gather all center locations
        minPair=min(loc,[],2)';
        maxPair=max(loc,[],2)';
        eset.executeExperimentFunction('trimTracks',[0 endSec+100],[minPair+trimBuffer maxPair-trimBuffer]);
        ao.trim = true;
    end

    if ~isfield(ao,'clean') || ao.clean==false
        ecl = ESetCleaner;
        ecl.minHTValid= .98;
        ecl.minDist=75;
        ecl.minSpeed = 0.75;
        ecl.minPts = 500;
        
        disp('cleaning eset');
        ecl.clean(eset);
        ao.clean=true;
    end
    
    %close figures opened by ecl.clean(eset)
    close all

    %write to file
    writeL=input('Do you want to write to mat file? 1 or 0: ');
    if writeL
        if~exist('toFile','var')|| exist('reSave','var') && reSave ||  exist('resetAll','var')
            disp('writing to mat file');
            
%             %gather all expt names
%             names=zeros(1,length(eset.expt));
%             for j=1:length(eset.expt)
%                 names(j)=eset.expt(j).fname;
%             end
%             
%             indFn=strfind(eset.expt(1).fname,'.');
%             fn=[eset.expt(1).fname(1:(indFn)-1) '_trim_' num2str(endSec) '_clean_' num2str(trimBuffer)];

            fname=eset.expt(1).fname;
            inds=strfind(fname, '\');
            fstub=[fname(1:inds(end)) 'MatFiles\'];
            mkdir(fstub);
            eset.toMatFiles([fstub 'DirectionalGradient_trim_' num2str(endSec) '_clean_' num2str(trimBuffer)]);
        end
    end
end   

% %fixes head/tail orientation
% 
% if ~isfield(ao,'fixht') || ao.fixht==false
%     disp('fixing head/tail orientation');
%     eset.executeTrackFunction('fixHTOrientation');
%     ao.fixht=true;
% end

%set headsweep theta min to 15 deg for segmentation (found this doesn't
%miss as many headsweeps for second instars
for j=1:length(eset)
    eset.expt(j).so.headswing_start=deg2rad(15);
end

% Set segmentation threshold speeds
if ~isfield(ao,'setSS') || ao.setSS==false
    disp('setting segment speeds');
    eset.executeTrackFunction('setSegmentSpeeds');
    ao.setSS=true;
end

% Segment the tracks
if ~isfield(ao,'seg') || ao.seg==false
    disp('segmenting tracks');
    eset.executeTrackFunction('segmentTrack');
    ao.seg=true;
end

%calculate conversion factor from frames to minutes
secPerFrame=eset(1).expt(1).dr.interpTime;
frame2min=(60/secPerFrame);
    
% Need to check segmentation by playing some movies here
if ~isfield(ao,'ao.checkSeg');
    checkSeg=input('Do you want to check segmentation? Enter 1 or 0: ');
    if checkSeg
        %use below for watching movies in order of ihtV score
        %figure;
        %eset.makeHistogram('ihtValid',0.8:0.01:1,'mean');
        %ihtV = eset.gatherField('ihtValid','mean');
        %[~,I] = sort(ihtV);
        %figure;
        for j=1:length(eset.expt)
            for i=1:length(eset.expt(j).track)
                eset.expt(j).track(i).playMovie;
                pause
            end
        end
    end
    ao.checkSeg=true;
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize fignum to zero
if checkSeg
    close all %close segmentation checking figures
end

fignum = 0;

%run marc's spatial navigation script
sno.angleBinSize=60;
sno.preferredDirection=180;
spatialNavigationMaggotFigures(eset,sno);

writeF=input('Do you want to write these figures to file?')
if writeF
    %get save directory
    fname=eset.expt(1).fname;
    inds=strfind(fname, '\');
    fstub=[fname(1:inds(end)) 'MatFiles\'];
    if exist(fstub,'file')~=7
        mkdir(fstub)
    end
    %print figures to file
    figHandles=findobj('Type','figure');
    count=1;
    for n=1:length(figHandles)
        figure(figHandles(n));
        print('-depsc',[fstub 'DirectionalGradient_figure' num2str(count)]);
        count=count+1;
    end
end



 