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

%plots path


%plots smoothed path


%1b) create a new figure, and plot the speed vs. time


%1c) make a histogram of the speed


% the full function call for the plotPath command is
% function plotPath(track, pathType, linetype, varargin)
%thus, to plot the smoothed location in blue, we would type
%track.plotPath('sloc', 'b-'); %this should look the same as your answer to 1a)

%1d) plot the displacement in red


%we can pass plotPath parameter/value pairs that can be passed to plot
%one such pair is 'LineWidth', lw;  which changes the width of the line we
%plot

%1e) plot the displacement in red, with a linewidth of 3



%% QUESTION 2



% segment the track into runs & reorientations

%these are the default segmentation options:
so = MaggotSegmentOptions;
%for an explanation of these parameters, type
%doc MaggotSegmentOptions;

%now, we want to set the segmentation speed on a case by case basis
%use the MaggotTrack command "setSegmentSpeeds" to do this
track.setSegmentSpeeds();

%2a) after you set segment speeds, find the new start & stop speeds


%2b) find the indices of all (interploated) points where the speed is less than 
%    the run stop speed

%another paremeter value pair you can pass to plotPath is 'highlightinds',
%inds;  this will cause little red dots to be drawn at every point marked by inds
%type help Track.plotPath for full details

%2c) plot the smoothed position of the track in blue, with all the points
%with speed < stop speed marked with red dots


%2d) to segment a single track, we just type
track.segmentTrack();

%after this runs, tell me how many runs, reorientations, and headSweeps are
%found in the track


%2e) close all open windows, then play a movie showing the segmentation,
%using track.playMovie
