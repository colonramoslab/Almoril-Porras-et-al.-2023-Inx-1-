function esetout = processLinjiaoFilesToMatfiles(dirname, varargin)
ts1 = tic;
minpts = 25;
frameDiff = 4; % stitch together tracks if first ended 2 or fewer frames before second started
maxDist = 0.1; %one mm
default_deltaT = 0.5;
varargin = assignApplicable(varargin);

[ps, genotype] = fileparts(dirname);
[~,condition] = fileparts(ps);

esetname = [genotype '_' condition];

if (~exist ('cc_setup1', 'var') || ~exist('cc_setup2', 'var'))
    load ('E:\worm thermotaxis bin files\calibration setup 1 20110309\camcal.mat');
    load ('E:\worm thermotaxis bin files\calibration setup 2 20110309\camcal.mat');
end

di = getDirectoryInformation(dirname);

setup2inds = [di.setupNumber] == 2;    
setup1inds = ~setup2inds;

cc = repmat(cc_setup1, size(di));
cc(setup2inds) = cc_setup2;
bfn = {di.binFileName};
eset = ExperimentSet.fromFiles(bfn{:}, 'camcalinfo', cc, 'minpts', minpts, 'parallel', false);
for j = 1:length(eset.expt)
    eset.expt(j).dr.interpTime = 0.25;
    tic
    [eset.expt(j).track.dr] = deal(eset.expt(j).dr);
    toc
end
eclnukespots = ESetCleaner();

eclnukespots.minDist = 0.1;
eclnukespots.askFirst = false;

eset.addTimingByFrameDifference(default_deltaT);    
eclnukespots.clean(eset);
eset.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);
ecl = ESetCleaner();
ecl.minDist = 1; %minimum distance 1 cm
ecl.minSpeed = 0.005; %minimum average speed 50 microns/sec
ecl.minPts = 100;

ecl.askFirst = false; 
ecl.clean(eset);

%{
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
eset.executeExperimentFunction('trimTracks', [], trimrect);
%}
disp('done with loading, stitching and cleaning');
toc(ts1)
disp ('basic calucations');
fields2calc = {'iloc', 'sloc','speed', 'scovRatio', 'covRatio', 'theta', 'deltatheta', 'displacement'};
for j = 1:length(fields2calc)
    try 
        eset.gatherField(fields2calc{j});
    catch me
        disp(me.getReport);
    end
end
        
toc(ts1);

mkdir (fullfile(dirname, 'matfiles'));
eset.toMatFiles(fullfile(dirname, 'matfiles',esetname));
disp('saved file');
toc(ts1)
if (nargout > 0)
    esetout = eset;
end
toc(ts1);