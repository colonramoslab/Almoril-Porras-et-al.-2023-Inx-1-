%%%%%%%%%%%%%%%%%%% BASIC LOCOMOTORY ANALYSIS VERSION 2 %%%%%%%%%%%%%%%%%%%
%
% Script loads data from a pair of .bin and .tim files; cleans tracks that
% are too short (in time or distance), too slow, or where the head/tail
% directions can't be determined very well; segments the tracks into runs
% and reorientations.  
%
% Then finds the mean larval speed, turning rate, head swing number, and 
% head swing angle.  Also prints the mean of these quantities for each
% track and generates a standard deviation and standard error accordingly. 
%
% [Last updated 08/16/2010]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mm / pixel conversion.  Use the values here or update depending on camera
% and plate position.  The easiest way to obtain this is to take a picture
% of a ruler sitting on the agar plate.  
lengthPerPixel = 1;
lengthPerPixel = 0.118; % (with the 8X lens)
%lengthPerPixel = 0.0766; % (with the 12X lens)
%
% Camera acquisition rate.  This is typically 4 frames/second.
framesPerSecond = 4;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PART I: LOAD AND PREPARE DATA %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load experiment set; don't load if one already has been loaded
if (~exist('eset','var'))
    % Loads data from the .bin file, and stores it in an ExperimentSet 
    % called "eset".  Asks the user for the location of the .bin file.  
    eset = ExperimentSet.fromFiles();
    
    % Defines an ESetCleaner called "ecl".  
    ecl = ESetCleaner;
    
    %use this command to get report when analyzing data for the first time:
    %ecl.getReport(eset);
end


% Clean up "bad" tracks
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    % Create an ESetCleaner called ecl.  (Redundant with the command above)
    ecl = ESetCleaner;
    
    % Properties of ecl.  These decide which tracks are thrown out.
    %
    ecl.minHTValid = .90; % percentage of track where motion is in the same
                          % direction as a tail-to-head vector
    ecl.minDist = 50;     % distance (in pixels) that a track must span
    ecl.minSpeed = 1;     % average speed cutoff, in pixels/s (not mm/s)
    ecl.minPts = 500;     % minimum number of frames in the track 
                          % (e.g., 500 frames = 125 seconds at 4 Hz acq.)
    
    % This command performs the actual track cleaning.  
    ecl.clean(eset);
    
    cleanEset = false;    % reset this to true if you want to re-clean    
end

% Fix head-tail orientation
existsAndDefault('fixht','true');
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

% Determine the speed thresholds used to decide when (1) a run ends and a
% new reorientation starts and (2) when a reorientation ends and a new run
% starts.  (These are not the same speeds).  These are found by looking at
% body orientation as a function of speed.  
existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% Segment the tracks (i.e., break them down into runs and reorientations)
existsAndDefault('segment', true);
if (segment)
    eset.executeTrackFunction('segmentTrack');
    segment = false; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PART II: ANALYZE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of tracks in the experiment (after the above cleaning)
numTracks = length(eset.expt.track);

% Define a track weighting vector, to show how much weight to give the
% result from each track.  Literally this is a vector that gives the
% duration (in frames) of each track.  Also generate a normalized version.
% Will use other weights later, depending on the quantity of interest.  
trackWeightFrame = zeros(1,numTracks);
for j=1:numTracks
   trackWeightFrame(j) = eset.expt.track(j).endFrame - eset.expt.track(j).startFrame + 1;    
end
trackWeightFrameNorm = trackWeightFrame/sum(trackWeightFrame); 

%%%%% A. SPEED %%%%%

% Find overall average speed.  Frames (track points) that occur during a
% reorientation event do not count when determining speed.  
sp = eset.gatherField('speed','run');
sp = lengthPerPixel * sp; % (convert speed to mm/s)
meanSpeed = mean(sp);
stdevSpeed = std(sp);

% Find the speed of each track separately.  
trackSpeed = zeros(1,numTracks);
for j=1:numTracks
    trackSpeed(j) = mean(eset.expt.track(j).getDerivedQuantity('speed',false,'run'));
end
trackSpeed = trackSpeed * lengthPerPixel; % (convert to mm/s)
wtSpeed = sum(trackWeightFrameNorm.*trackSpeed);

% Find weighted standard deviation and standard error of larval speed.  
wtSDspeed = sqrt( sum( trackWeightFrameNorm.*(trackSpeed-wtmeanSpeed).*(trackSpeed-wtmeanSpeed) ) );
wtSEspeed = wtSDspeed / sqrt(numTracks-1);


%%%%% B. TURN RATE %%%%%

% Find the turning rate for each track.  
trackTurnrate = zeros(1,numTracks);
for j=1:numTracks
    % For each track, count the number of turns (reorientations with at
    % least one head sweep).  
    numTurns=0;
    for k=1:length(eset.expt.track(j).reorientation)
        if (eset.expt.track(j).reorientation(k).numHS>0)
            numTurns = numTurns + 1;
        end
    end
    
    % Assign the number of turns to the turn rate vector.  Obtain a turn
    % rate by dividing by the number of points in the track.  
    trackTurnrate(j) = numTurns/(eset.expt.track(j).npts);
end

% Convert turns/frame into turns/minute
trackTurnrate = trackTurnrate * framesPerSecond * 60;

% Find average reorientation rate per larva, weighted by track duration.  
wtTurnrate = sum(trackWeightFrameNorm .* trackTurnrate);

% Find weighted standard deviation and standard error of larval turn rate.  
wtSDturnrate = sqrt( sum( trackWeightFrameNorm.*(trackTurnrate-wtTurnrate).*(trackTurnrate-wtTurnrate) ) );
wtSEturnrate = wtSDturnrate / sqrt(numTracks-1);


%%%%% C. NUMBER OF HEAD SWEEPS %%%%%

% Overall average number of head sweeps (per reorientation)
numberHS = eset.gatherSubField('reorientation','numHS');
meanHS = mean(numberHS(numberHS>0));
stdevHS = std(numberHS(numberHS>0));

% Find the number of head sweeps for each track.  Also collect the number 
% of reorientations in each track (for weighting).
trackNumReo = zeros(1,numTracks);
trackNumHS = zeros(1,numTracks);
for j=1:numTracks
   
    for k=1:length(eset.expt.track(j).reorientation)
        % Number of head sweeps in the kth reorientation of the jth track
        numHSjk = eset.expt.track(j).reorientation(k).numHS;
        
        % If the reorientation has head sweeps, at that number to the
        % total, and also add to the reorientation count.  
        if(numHSjk>0)
            trackNumReo(j) = trackNumReo(j) + 1;
            trackNumHS(j) = trackNumHS(j) + numHSjk;
        end
    end   
end

% New weighting scheme based on the number of reorientations in the track,
% instead of on the number of frames, as was done above.  
trackWeightReo = trackNumReo;
trackWeightReoNorm = trackWeightReo/sum(trackWeightReo);

% Have the total number of head sweeps in each track.  Convert to head
% sweeps per reorientation in each track:
trackNumHSreo = trackNumHS ./ trackNumReo;

% Find average number of head sweeps per reorientation, weighted by the
% number of reorientations in the track.  (This should exactly equal the
% meanHS calculated above, unlike that for speed or reorientation rate).  
wtNumHS = sum(trackWeightReoNorm .* trackNumHSreo);

% Find weighted standard deviation and standard error of head sweep number.  
wtSDnumHS = sqrt( sum( trackWeightReoNorm.*(trackNumHSreo-wtNumHS).*(trackNumHSreo-wtNumHS) ) );
wtSEnumHS = wtSDnumHS / sqrt(numTracks-1);


%%%%% D. HEAD SWEEP DEPTH %%%%%

% Overall average depth of a head sweep
sweepDepth = eset.gatherSubField('headSwing', 'maxTheta');
sweepDepth = abs(sweepDepth)*(180/pi);
meanSweepDepth = mean(sweepDepth);
stdevSweepDepth = std(sweepDepth);

% Find the average head sweep depth for each track.  Also, for fun, find
% the average depth of just the first head sweep for a given reorientation.
trackHSdepth = zeros(1,numTracks);
trackFirstHSdepth = zeros(1,numTracks);
for j=1:numTracks
    
    for k=1:length(eset.expt.track(j).reorientation)
        % Number of head sweeps in the kth reorientation of the jth track
        numHSjk = eset.expt.track(j).reorientation(k).numHS;
        
        if(numHSjk>0)
            
            for m=1:numHSjk
                
                thetajkm = abs(eset.expt.track(j).reorientation(k).headSwing(m).maxTheta)*(180/pi);
                trackHSdepth(j) = trackHSdepth(j) + thetajkm;
                
                if(m==1)
                    trackFirstHSdepth(j) = trackFirstHSdepth(j) + thetajkm;
                end
            end
        end
    end
end

% Another weighting scheme, based on the number of head sweeps in each
% track.  Results from this should look very similar to those with weights
% based on reorientations. 
trackWeightHS = trackNumHS;
trackWeightHSnorm = trackWeightHS/sum(trackWeightHS);

% Have the sum of all the head sweep angles.  Convert to angle per head
% swing.  
trackHSdepth = trackHSdepth ./ trackNumHS;
trackFirstHSdepth = trackFirstHSdepth ./ trackNumReo;  

% Find the average angle of a head sweep, weighted by the number of
% head sweeps in the track.    
wtSweepdepth = sum(trackWeightHSnorm .* trackHSdepth);

% Find weighted standard deviation and standard error of head sweep depth.  
wtSDsweepdepth = sqrt( sum( trackWeightHSnorm.*(trackHSdepth-wtSweepdepth).*(trackHSdepth-wtSweepdepth) ) );
wtSEsweepdepth = wtSDsweepdepth / sqrt(numTracks-1);

% Average angle of the first head sweep in a reorientation, weighted by the
% number of reorientations in teh track.  Also find standard deviation and
% standard error.  
wtFirstSweepdepth = sum(trackWeightReoNorm .* trackFirstHSdepth);
wtSDfirstSweepdepth = sqrt( sum( trackWeightReoNorm.*(trackFirstHSdepth-wtFirstSweepdepth).*(trackFirstHSdepth-wtFirstSweepdepth) ) );
wtSEfirstSweepdepth = wtSDfirstSweepdepth / sqrt(numTracks-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PART III: DISPLAY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('***************************************************************');

display(['Analyzed ' num2str(numTracks) ' tracks']);

display('***************************************************************');

display('SPEED');

display(['  average larval speed: ' num2str(wtSpeed,3) ' mm/s']);
display(['  (standard deviation ' num2str(wtSDspeed,3) ')']);
display(['  (standard error ' num2str(wtSEspeed,3) ')']);

display('***************************************************************');

display('TURN RATE');

display(['  average turning rate per larva: ' num2str(wtTurnrate,3) ' turns/min']);
display(['  (standard deviation ' num2str(wtSDturnrate,3) ')']);
display(['  (standard error ' num2str(wtSEturnrate,3) ')']);

display('***************************************************************');

display('HEAD SWEEP NUMBER');

display(['  average number of head sweeps per larval reorientation: ' num2str(wtNumHS,3)]);
display(['  (standard deviation ' num2str(wtSDnumHS,3) ')']);
display(['  (standard error ' num2str(wtSEnumHS,3) ')']);

display('***************************************************************');

display('HEAD SWEEP ANGLE');

display(['  average head sweep angle per larval head sweep: ' num2str(wtSweepdepth,3)]);
display(['  (standard deviation ' num2str(wtSDsweepdepth,3) ')']);
display(['  (standard error ' num2str(wtSEsweepdepth,3) ')']);

display(['  average first head sweep angle per larval head sweep: ' num2str(wtFirstSweepdepth,3)]);
display(['  (standard deviation ' num2str(wtSDfirstSweepdepth,3) ')']);
display(['  (standard error ' num2str(wtSEfirstSweepdepth,3) ')']);

display('***************************************************************');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






