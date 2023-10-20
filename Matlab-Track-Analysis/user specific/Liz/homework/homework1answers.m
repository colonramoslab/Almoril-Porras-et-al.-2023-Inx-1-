%% HOMEWORK #1 
% WORKING WITH TRACKS

%% SETTING UP
% preliminary stuff to get you ready to go

%first, let's load an experiment from disk;
fn = '\\labnas1\Share\Phototaxis\Data Extracted\Gradient\VaryGradient\26_1515_tracks.bin'; %relpace with the name of your favorite bin file

if (~exist('eset','var')) % only load if you haven't already
    eset = ExperimentSet.fromFiles(fn);
end

%now let's pick out a track to work with;  (pick a nice long one)
%plot some tracks to pick a long one

tnum = 2; %you can change this
track = eset.expt(1).track(tnum);

%%
%QUESTION 1
% plot the track path & parameters
%1a) create a new figure, and plot the smoothed path of the track
fignum = 0;
%plots path
fignum = fignum + 1; figure(fignum); clf(fignum);
track.plotPath;
title('path');
axis equal;

%plots smoothed path
fignum = fignum + 1; figure(fignum); clf(fignum);
sloc=track.getDerivedQuantity('sloc');
plot(sloc(1,:), sloc(2,:));
title('smoothed path');
axis equal;

%1b) create a new figure, and plot the speed vs. time
fignum = fignum + 1; figure(fignum); clf(fignum);
sp=track.getDerivedQuantity('speed');
t=track.getDerivedQuantity('eti');
plot(t,sp);
title('speed vs. time');
xlabel('time');
ylabel('speed');

%1c) make a histogram of the speed
fignum = fignum + 1; figure(fignum); clf(fignum);
hist(sp);
title('speed hist');

% the full function call for the plotPath command is
% function plotPath(track, pathType, linetype, varargin)
%thus, to plot the smoothed location in blue, we would type
%track.plotPath('sloc', 'b-'); %this should look the same as your answer to 1a)

%1d) plot the displacement in red
fignum = fignum + 1; figure(fignum); clf(fignum);
track.plotPath('displacement','r-');

%we can pass plotPath parameter/value pairs that can be passed to plot
%one such pair is 'LineWidth', lw;  which changes the width of the line we
%plot

%1e) plot the displacement in red, with a linewidth of 3
fignum = fignum + 1; figure(fignum); clf(fignum);
track.plotPath('displacement','r-','LineWidth',3);



%% QUESTION 2

fignum=0;

% segment the track into runs & reorientations

%these are the default segmentation options:
so = MaggotSegmentOptions;
%for an explanation of these parameters, type
%doc MaggotSegmentOptions;

%now, we want to set the segmentation speed on a case by case basis
%use the MaggotTrack command "setSegmentSpeeds" to do this
track.setSegmentSpeeds();

%2a) after you set segment speeds, find the new start & stop speeds
disp('start speed cut:'); disp(track.so.start_speed_cut);
disp('stop speed cut:'); disp(track.so.stop_speed_cut);

%2b) find the indices of all (interploated) points where the speed is less than 
%    the run stop speed
s = track.getDerivedQuantity('speed');
inds = find( s < track.so.stop_speed_cut);

%another paremeter value pair you can pass to plotPath is 'highlightinds',
%inds;  this will cause little red dots to be drawn at every point marked by inds
%type help Track.plotPath for full details

%2c) plot the smoothed position of the track in blue, with all the points
%with speed < stop speed marked with red dots
fignum = fignum + 1; figure(fignum); clf(fignum);
track.plotPath('sloc', 'b-', 'highlightinds', inds);

%2d) to segment a single track, we just type
track.segmentTrack();

%after this runs, tell me how many runs, reorientations, and headSweeps are
%found in the track
disp('num runs:'); disp(length(track.run));
disp('num reorientations:'); disp(length(track.reorientation));
disp('num headsweeps:'); disp(length(track.headSwing));

%2e) close all open windows, then play a movie showing the segmentation,
%using track.playMovie
close all
track.playMovie;
