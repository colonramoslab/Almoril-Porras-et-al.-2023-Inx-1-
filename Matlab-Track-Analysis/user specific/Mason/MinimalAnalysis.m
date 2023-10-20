% This first part is copied from Liz's example script.  

% minimum track length (in *frames*) that will be included.  
%mintracklength=500;

% number of animals (or should we use number of tracks for standard error?)
% numMaggots = 30

% mm / pixel conversion
lengthPerPixel = 1;
lengthPerPixel = 0.118; % (with the 8X lens)
%lengthPerPixel = 0.0766; % (with the 12X lens)

%load experiment set, dont load if one already has been loaded
if (~exist('eset','var'))
    %eset=ExperimentSet.fromFiles('J:\VANVACTORLAB\06072010\Analysis\Movie6','minpts',mintracklength);
    % this part loads all the data from your .bin file, now stored in an
    % ExperimentSet called "eset".  
    eset=ExperimentSet.fromFiles();
    %use the below   code to get report when analyzing data for the first
    %time 
    % Defines an ESetCleaner called "ecl".  
    ecl = ESetCleaner;
    %ecl.getReport(eset);
end


%clean up bad tracks and treshold data
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid=.90;
    ecl.minDist=50;
    ecl.minSpeed = 1;% in pixels/s, not mm/s
    ecl.minPts = 500;
    
    ecl.clean(eset);
    
    %trimRect = [355 800 2145 1795];
    %eset.executeExperimentFunction('trimTracks', [], trimRect);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% speed
sp = eset.gatherField('speed','run');
sp = lengthPerPixel * sp; % converts to mm/s
meansp = mean(sp)
stdevsp = std(sp)

% turning rate...
runendtime = eset.gatherField('eti', 'runend');
eti = eset.gatherField('eti','run');
numRuns = length(runendtime); 
numFrames = length(eti);
turnrate = (numRuns/numFrames)*240 %240 is 4 frames/sec x 60 sec/min

% that's the overall turning rate.  also might be useful to find the
% turning rate for each track, then average those.  the overall mean will
% be the same, but this will give a better idea of the spread and allow us
% to put error bars on things.  
for j=1:length(eset.expt.track)
    numTurns=0;
    for k=1:length(eset.expt.track(j).reorientation)
        if (eset.expt.track(j).reorientation(k).numHS>0)
            numTurns = numTurns + 1;
        end
    end
    trackturnrate(j)=numTurns/(eset.expt.track(j).npts);
end
trackturnrate = trackturnrate * 240;
% This turn rate doesn't count reorientations with 0 head sweeps.  
% (this is the turn rate that was included in the spreadsheet).
turnrate2 = mean(trackturnrate)
stdevTurnRate = std(trackturnrate)

% avg. # head sweeps...
numberHS = eset.gatherSubField('reorientation','numHS');
meanHS = mean(numberHS(numberHS>0))
stdevHS = std(numberHS(numberHS>0))

% avg. depth of sweep
sweepDepth = eset.gatherSubField('headSwing', 'maxTheta');
sweepDepth = abs(sweepDepth)*(180/pi);
meanSweepDepth = mean(sweepDepth)
stdevSweepDepth = std(sweepDepth)

% number of tracks
numTracks = length(eset.expt.track)



