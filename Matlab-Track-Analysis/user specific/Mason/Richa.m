% This first part is copied from Liz's example script.  

%first example, loading a single experiment
mintracklength = 500;

% number of animals (or should we use number of tracks for standard error?)
% numMaggots = 30

% mm / pixel conversion
lengthPerPixel = 0.118 % (with the 8X lens)
lengthPerPixel = 0.0766 % (with the 12X lens)

%load experiment set, dont load if one already has been loaded
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles('C:\Documents and Settings\Mason\Desktop\RichaData\05282010\movie4','minpts',mintracklength);
    %use the below code to get report when analyzing data for the first
    %time
    ecl = ESetCleaner;
    %ecl.getReport(eset);
end


%clean up bad tracks and treshold data
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid=.90;
    ecl.minDist=0;
    ecl.minSpeed = 1;
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

% speed
sp = eset.gatherField('speed','run');
sp = lengthPerPixel * sp
meansp = mean(sp)
stdevsp = std(sp)
numRuns = length(runendtime)

% turning rate...
runendtime = eset.gatherField('eti', 'runend');
eti = eset.gatherField('eti','run');
turnrate = (length(runendtime)/length(eti))*240

% avg. # head sweeps...
numberHS = eset.gatherSubField('reorientation','numHS');
meanHS = mean(numberHS)
stdevHS = std(numberHS)

% avg. depth of sweep
sweepDepth = eset.gatherSubField('headSwing', 'maxTheta');
sweepDepth = abs(sweepDepth);
meanSweepDepth = mean(sweepDepth)
stdevSweepDepth = std(sweepDepth)


