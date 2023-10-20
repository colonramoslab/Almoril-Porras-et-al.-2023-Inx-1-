%% ANALYZING WORM FILES
% A Step by Step guide to using the track analysis code
%%

% This is a starter list of things to do with worm experiment files
% It was created for Linjiao, but feel free to add other commands with 
% documentation

%% LOADING FILES FROM DISK
% how to load in your experiments

% To load experiment sets from files, use the ExperimentSet class
% for information on the ExperimentSet class, type
doc ExperimentSet

% To load all tracks at least 200 points long from a single file
eset = ExperimentSet.fromFiles(filename,'minpts', 200);
% if you do not specify a minimum length, the default is 100 points

% To load tracks from multiple files
eset = ExperimentSet.fromFiles(filename1, filename2, filename3, 'minpts', 200);

% To load tracks from all .bin files in a directory
eset = ExperimentSet.fromFiles(directoryname, 'minpts', 200);

% To select tracks by hand with GUI
eset = ExperimentSet.fromFiles()
% or
eset = ExperimentSet.fromFiles('minpts', 200);

% All the information from a single file is stored in an Experiment object
% For information on experiments, type
doc Experiment

% To execute the same function for all experiments in an experiment set, use
eset.executeExperimentFunction('functionName', arglist);

%% CLEANING UP FILES
% Optional arguments to stitch together tracks, remove spurious points

%%% STITCHING TRACKS
% If the track extraction software missed the animal in one frame, you'll
% have two tracks instead of one;  here's how to stitch them back together

% To stitch together tracks in an experiment (to clean up missing frames)
doc Experiment.stitchTracks
eset.expt(j).stitchTracks(maxFrameInterval, maxDist)
% to stitch all tracks in an experiment set
eset.executeExperimentFunction('stitchTracks', maxFrameInterval, maxDist);

% each time the program finds tracks to stitch, it will show you a graph (in
% figure 1) of the two tracks, with the starting track in green, the
% track to be merged onto the end in red, and any other nearby tracks in gray
% in figure 2, it will play a movie of the end of track 1 followed by the
% start of track 2
% you will get a prompt asking stitch tracks y/n
% to stitch the tracks, type 'y' and hit enter;  to leave them separate 'n'
% and enter.  to see the movie again, just hit enter

% to turn this feature off and have all track matches stitched together
% automatically, 
eset.expt(j).stitchTracks(maxFrameInterval, maxDist, 'interactive', false);
eset.executeExperimentFunction('stitchTracks', maxFrameInterval, maxDist, 'interactive', false);

%%% CLEANING UP BAD TRACKS
% To clean up short tracks, tracks where the worm just sits there, etc.

% for information on the ESetCleaner class
doc ESetCleaner

ecl = ESetCleaner();
ecl.minSpeed = xx; %minimum average speed of the track to be kept
ecl.minDist = xx; %minimum distance must travel from the starting point to be kept
ecl.minPts = xx; %minimum number of points to be kept

% generate graphs showing the distributions of speed, dist, and pts and
% which would be cut
ecl.getReport (eset);

% actually do the trimming (will show you report and ask if you want to
% proceed)
ecl.clean(eset);

%%% REMOVING EXCESS TRACK
% To get rid of portions of track that leave a certain area (e.g. to ignore
% worms that reach the edge & then return), see help file:
doc Experiment.trimTracks

% To remove tracks that start outside a valid time or space range, see help
% file:
doc Experiment.pruneTracks

% To execute the trim or prune on all experiments in the experiment set
% using the same parameters
eset.executeExperimentFunction('trimTracks', timerange, validrect)
eset.executeExperimentFunction ('pruneTracks', starttimerange, startrect)
 
%% SEGMENTING TRACKS INTO RUNS AND REORIENTATIONS
% how to parse the track into periods of runs and reorientation

%%% SEGMENTATION OPTIONS
% These set the rules for segmenting tracks

% from the help information:
doc WormSegmentationOptions

% first sharp turns are flagged as any point at which the velocity
% direction is changing faster than DTHETATHRESH
% any sharp turns that are less than JOINSTPTS apart are merged into a
% single sharp turn
% the velocity angle into a sharp turn is the angle PTBUFFER points before
% the turn; and the velocity angle out is PTBUFFER points after
% if the in and out velocity angles are within ALIGNEDTHETA, the sharp
% turn is flagged as a reversal (if out is ~180 degrees from in) or a
% blip (if in the same direction)
% otherwise, flagged as an omega turn
% all sharp turns within MINRUNTIME (in seconds) of each other are
% grouped into a single reorientation
% a reorientation ends at the point where dtheta/dt is less than
% STRAIGHTTHETATHRESH after the last sharp turn in that reorientation
 
%%% SEGMENTING TRACKS
% Doing the segmentation

% to use the same set of rules for all tracks
eset.expt(j).segmentTracks(worm_segment_options)
% or
eset.executeExperimentFunction('segmentTracks', worm_segment_options);

% to use each a different segment option for each track, modify the track's
% segment options individually,
eset.expt(j).track(k).so = new_segment_options;
% then
eset.expt(j).segmentTracks()
% or
eset.executeExperimentFunction('segmentTracks');

%% SPECIFIC EXAMPLES
% These are commands to do specific things Linjiao requested

%% plot of all the tracks with the start of tracks superposed at the center
% use the 'displacement' field

% two possibilities:
% method 1:
figure(fignum); clf(fignum);  %create a figure and clear it
hold on;  
eset.executeTrackFunction('plotPath', 'displacement'); % plot all paths
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
hold off;

% method 2:
figure(fignum); clf(fignum);  %create a figure and clear it
alltracks = [eset.expt.track]; %get the handles for all tracks in one vector
alltracks.plotPath('displacement'); % plot all paths
hold on
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
hold off

%% x position of the center of the mass over time
% use the 'sloc' (smoothed location) field, then subtract off start point

timeAxis = firstTime:DT:lastTime; % in seconds
[time, centerofmass] = eset.meanField2vsField1('eti', 'sloc', timeAxis);
xposition = centerofmass(1,:);
plot (time, xposition - xposition(1,1));

%% x direction drift velocity over time
% this is just the derivative of the center of mass over time

% run previous commands to find center of mass, then
sigma = 1; %amount of smoothing
xcomvel = deriv(xposition)./deriv(time);
plot (time, xcomvel);

%% scattering plot of run duration at different directions
% get the two fields you're interested in and plot them against each other

rundirection = eset.gatherSubField('run', 'meanTheta');
runduration = eset.gatherSubField('run', 'runTime');
plot (rad2deg(rundirection), runduration);

%% histogram of run direction
% use the makeSubFieldHistogram command

% generate a theta axis; although run angles run between -pi and pi
% radians, using the 'polar' option below shifts them to between 0 and 2 pi
thetaAxis = deg2rad(0:30:330);
eset.makeSubFieldHistogram('run', 'meanTheta', thetaAxis, 'polar', true, 'r2d', true);


%% run speed at different directions
% use the meanField2vsField1 commanda again

% here we need the theta axis to run between -pi and pi
thetaAxis = deg2rad(-180:30:180);
[newtheta, meanspeed] = eset.meanField2vsField1('theta', 'speed', thetaAxis, 'run');
plot (rad2deg(newtheta), meanspeed);

%% excluding endpoints from runs
thetaAxis = deg2rad(-180:30:180);
%to exclude first 5 and last 5 pts from run
runtheta = eset.gatherFromSubField('run', 'theta', 'trimpts', 5);
runspeed = eset.gatherFromSubField('run', 'speed', 'trimpts', 5);
[newtheta, meanspeed, stderrorspeed] = meanyvsx(runtheta, runspeed, thetaAxis);
errorbar(rad2deg(newtheta), meanspeed, stderrorspeed);


%% scattering plot of run direction change at different initial run direction
% gather the start and end directions then process

runstart = eset.gatherSubField('run', 'startTheta');
runend = eset.gatherSubField('run', 'endTheta');
rundt = diff(unwrap([runstart;runend])); % this takes care of 2 pi discontinuities

% group into quadrants
startedRight = abs(runstart) < pi / 4;  % etc.

% make histograms by quadrant
dthetaAxis = deg2rad(-180:15:180);
plot (rad2deg(dthetaAxis), hist(rundt(startedRight), dthetaAxis)); % etc.
