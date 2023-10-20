function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%loads, stitches, trims, cleans, etc. then saves to mat files in a
%'matfiles' subdirectory 
%varargin: minpts, frameDiff, buffer, trimrect, default_deltaT,
%checker_calc, timerange

minpts = 50;
frameDiff = 2; % stitch together tracks if first ended 2 or fewer frames before second started
if (isempty(camcalinfo))
    maxDist = 10; % stitch together tracks if first ended within 10 pixels of second's start
else
    maxDist = 0.1; %one mm
end
buffer = [];
trimrect = [];
default_deltaT = 0.2;

checker_calc = [];
timerange = [];
varargin = assignApplicable(varargin);
eclexisted = existsAndDefault('ecl', ESetCleaner());
existsAndDefault('camcalinfo', []);

eclnukespots = ESetCleaner();
eclnukespots.minHTValid = 0.5;
if (~isempty(camcalinfo)) 
    eclnukespots.minDist = 0.1;
    eclnukespots.minHTValid = 0.5;
    
end
eclnukespots.askFirst = false;
ts1 = tic;

eset = ExperimentSet.fromFiles(basedir, 'minpts', minpts, 'camcalinfo', camcalinfo, 'parallel', true, varargin{:});
eset.addTimingByFrameDifference(default_deltaT);
    
eclnukespots.clean(eset);
if (~eclexisted ) 
    if (~isempty(camcalinfo))
    %real units, assume cm
        ecl.minDist = 0.1; %minimum distance 1 mm
        ecl.minSpeed = 0.01; %minimum average speed 100 microns/sec
    else
        ecl.minDist = 10; %pixels
        ecl.minSpeed = 0.1; %pixels/second
    end
    ecl.minPts = ceil(30 / eset.expt(1).dr.interpTime);
end

eset.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);
ecl.askFirst = false; 
ecl.showFigsInReport = false;
ecl.getReport(eset);
ecl.clean(eset);

if (isempty(trimrect)) 
    il = eset.gatherField('iloc');
    if (isempty(buffer))
        if(isempty(camcalinfo))
            buffer = 25; %pixels
        else
            buffer = 0.15; %cm
        end
    end
    ll = min(il,[],2) + buffer;
    ur = max(il,[],2) - buffer;
    trimrect = [ll(1) ll(2) ur(1) ur(2)];
end
eset.executeExperimentFunction('trimTracks', timerange, trimrect);
if (isa (eset.expt(1).track(1), 'MaggotTrack'))
    eset.executeTrackFunction('setSegmentSpeeds');
    eset.gatherField('sspineTheta');
    eset.gatherField('vel_dp');
end
disp('done with loading, stitching and cleaning');
toc(ts1)

if (~isempty(checker_calc))
    try
        disp ('assigning checker track data');
        checker_calc.assignGlobals(eset.expt);
    catch me
        disp(me.getReport());
    end
    toc(ts1);
end
ts1 = tic;
mkdir (fullfile(basedir, 'matfiles'));
eset.toMatFiles(fullfile(basedir, 'matfiles',esetname));
disp('saved file');
toc(ts1)

