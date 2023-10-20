%% Set up cryophilic point labeling
% An example script with annotations
%%

existsAndDefault('fromScratch', false);
%% LOADING FILES FROM DISK
% loading specific files by name
ts1 = tic;
basedir = 'E:\larvalco2\Extracted\50 mL CO2 in 2 L air\gr63a';
esetname = 'co2_gr63a';
if (fromScratch)
    close all;
      minpts = 50;

    % this code snippet loads the files if we haven't already loaded them, but
    % otherwise skips them; that way we can change the script and rerun it
    % without having to reload the files
    if (~exist('co2_gr63a', 'var'))
        co2_gr63a = ExperimentSet.fromFiles(basedir, 'minpts', minpts);
    end
    
    

    %% STITCH TRACKS
    % sometimes we miss a frame, so let's stitch together tracks that are close
    % by

    frameDiff = 2; % stitch together tracks if first ended 3 or fewer frames before second started
    maxDist = 10; % stitch together tracks if first ended within 7 pixels of second's start

    % For the script, I am executing this function with interactive off, but if
    % you set interactive to true, it will show you each potential stitch and
    % let you decide whether or not to stitch it
    %co2_gr63a.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);

    %% CLEAN UP TRACKS
    % get rid of any tracks that don't go anywhere

    % create an EsetCleaner object

    ecl = ESetCleaner();

    % now let's look at the autogenerated report
    % let's get rid of all tracks less than 750 points and speed less than 0.4
    % pixels per second
    ecl.minPts = 500;
    ecl.minSpeed = 2;
    ecl.minHTValid = 0.95;

    ecl.getReport(co2_gr63a);
    


    % we've already shown the report, so we don't need to have it ask us first,
    % for the purposes of this script;  generally a good idea to leave this
    % enabled
    ecl.askFirst = false; 

    ecl.clean(co2_gr63a);
    
    buffer = 25;
    trimrect = [buffer buffer 2492-buffer 1944-buffer];
    co2_gr63a.executeExperimentFunction('trimTracks', [], trimrect);
    
    
    disp('done with loading, stitching and cleaning');
    toc(ts1)
    ts1 = tic;
   % save (fullfile(basedir, 'marcmatfile.mat'), 'cryo');
    mkdir (fullfile(basedir, 'matfiles'));
    co2_gr63a.defaultTitle = 'Gr63a';
    co2_gr63a.toMatFiles(fullfile(basedir, 'matfiles',esetname));
    disp('saved file');
    toc(ts1)
    fromScratch = false;
    resegment = true;
else
    if (~exist('co2_gr63a', 'var'))
        ts1 = tic;
     
        %load (fullfile(basedir, 'marcmatfile.mat'));
        co2_gr63a = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles',esetname));
        disp ('loaded eset from mat file');
        toc(ts1)
        resegment = true;
    end
end

existsAndDefault('resegment', true);
if (resegment)
    co2_gr63a.executeTrackFunction('setSegmentSpeeds');
    co2_gr63a.executeTrackFunction('segmentTrack');
    resegment = false;
end
sno.angleBinSize = 60;
sno.preferredDirection = 180;
if (exist ('co2_4pct', 'var') && exist('co2_control', 'var'))
    ad = spatialNavigationMaggotAnalysis([co2_4pct, co2_control, co2_gr63a], sno);
    po(1).legendEntry = 'CO$_2$ 0-4$\%$';
    po(2).legendEntry = 'control';
    po(3).legendEntry = 'Gr63a';
    spatialNavigationMaggotFigures(ad, sno, po);
end
