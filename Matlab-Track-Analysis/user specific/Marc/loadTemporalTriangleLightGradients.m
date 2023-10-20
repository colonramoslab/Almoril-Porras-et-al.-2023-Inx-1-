existsAndDefault('fromScratch', false);
%% LOADING FILES FROM DISK
% loading specific files by name
ts1 = tic;
basedir = '\\LABNAS2\Phototaxis\Extracted Phototaxis Data\ReExtracted\Temporal Gradients\Triangle\';
times = {'100s', '200s', '400s', '800s'};
esetbase = 'triangle_';
if (fromScratch)
    for j = 1:length(times)
        close all;
        minpts = 50;
        esetname = [esetbase times{j}];
        % this code snippet loads the files if we haven't already loaded them, but
        % otherwise skips them; that way we can change the script and rerun it
        % without having to reload the files
        if (~exist(esetname, 'var'))
            esettemp = ExperimentSet.fromFiles(fullfile(basedir,times{j}), 'minpts', minpts);
        else
            continue;
        end
    
    

        %% STITCH TRACKS
        % sometimes we miss a frame, so let's stitch together tracks that are close
        % by

        frameDiff = 2; % stitch together tracks if first ended 3 or fewer frames before second started
        maxDist = 10; % stitch together tracks if first ended within 7 pixels of second's start

        % For the script, I am executing this function with interactive off,
        % but if
        % you set interactive to true, it will show you each potential stitch and
        % let you decide whether or not to stitch it
        esettemp.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);

        %% CLEAN UP TRACKS
        % get rid of any tracks that don't go anywhere

        % create an EsetCleaner object

        ecl = ESetCleaner();
        ecl.minPts = 500;
        ecl.minSpeed = 1.0;
        ecl.minHTValid = 0.95;
        ecl.askFirst = false; 

        ecl.clean(esettemp);
    
        allloc = esettemp.gatherField('iloc');
        minl = min(allloc, [], 2);
        maxl = max(allloc, [], 2);
    
        buffer = 25;
        trimrect = [minl(1)+buffer minl(2)+buffer maxl(1)-buffer maxl(1)+buffer];
        esettemp.executeExperimentFunction('trimTracks', [], trimrect);
    
    
        disp('done with loading, stitching and cleaning');
        toc(ts1)
        ts1 = tic;
   % save (fullfile(basedir, 'marcmatfile.mat'), 'cryo');
        mkdir (fullfile(basedir, times{j}, 'matfiles'));
        esettemp.toMatFiles(fullfile(basedir, times{j}, 'matfiles',esetname));
        disp('saved file');
        toc(ts1)       
        disp (['changing esettemp name to ' esetname]);
        ts1 = tic;
        eval([esetname ' = esettemp;']);
        clear esttemp;
        toc(ts1);
        memory
        fromScratch = false;
    end
else
    for j = 1:length(times)
        close all;
        minpts = 50;
        esetname = [esetbase times{j}];

        if (~exist(esetname, 'var'))
            ts1 = tic;
     
            %load (fullfile(basedir, 'marcmatfile.mat'));
            eval([esetname ' = ExperimentSet.fromMatFiles(fullfile(basedir, times{j}, ''matfiles'',esetname));']);
            disp (['loaded ' esetname ' from mat file']);
            toc(ts1)
        end
    end
end
if (~exist('alltriangles', 'var'))
    command = 'alltriangles = [';
    for j = 1:length(times)
        command = [command ' ' esetbase times{j}];
    end
    command = [command '];'];

    eval(command);
end

existsAndDefault('resegment', true);
if (resegment)
    for j = 1:length(alltriangles)
        alltriangles(j).executeTrackFunction('setSegmentSpeeds');
        alltriangles(j).executeTrackFunction('segmentTrack');
        resegment = false;
    end
end

existsAndDefault('addGlobalQuantities', true);
if (addGlobalQuantities)
    for j = 1:length(alltriangles)
        for k = 1:length(alltriangles(j).expt)
            tic
            gq = alltriangles(j).expt(k).globalQuantity(2);
            dlt = gq.yData;
            lt = alltriangles(j).expt(k).globalQuantity(1).yData;
            cycleStartLow = find(diff(dlt > 0) > 0);
            cycleStartHigh = find(diff(dlt < 0) > 0);
            eti = gq.xData;
            eti(1:(cycleStartLow(1) - 1)) = -100;

            for m = 1:length(cycleStartLow)
                eti(cycleStartLow(m):end) = eti(cycleStartLow(m):end) - eti(cycleStartLow(m));
            end
            gq.fieldname = 'cycleTimeLow';
            gq.yData = eti;
            alltriangles(j).expt(k).addGlobalQuantity(gq);

            eti = gq.xData;
            eti(1:(cycleStartHigh(1) - 1)) = -100;

            for m = 1:length(cycleStartHigh)
                eti(cycleStartHigh(m):end) = eti(cycleStartHigh(m):end) - eti(cycleStartHigh(m));
            end
            gq.fieldname = 'cycleTimeHigh';
            gq.yData = eti;
            alltriangles(j).expt(k).addGlobalQuantity(gq);
            
            gq.fieldname = 'dloglt';
            ltmin = min(lt(lt > 0));
            lt (lt <= 0) = ltmin;
            gq.yData = dlt./lt;
            alltriangles(j).expt(k).addGlobalQuantity(gq);
            toc
        end
    end
    addGlobalQuantities = false;
end
        