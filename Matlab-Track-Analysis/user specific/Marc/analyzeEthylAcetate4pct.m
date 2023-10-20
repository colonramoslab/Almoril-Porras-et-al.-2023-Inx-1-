%this is the directory where your .bin and .tim files are stored
basedir = 'd:\marc Processed\maggots\ethyl acetate 4 pct 20 2000\';

minLoadPts = 50; %don't load any tracks shorter than this

%load all the files, this will take a while, so we only do it
%if we haven't already loaded them
if (~exist('odor4pct', 'var'))
    odor4pct = ExperimentSet.fromFiles(basedir,'minpts',minLoadPts); 
    disp('files loaded; fixing head tail orientation; this may take a while');
    tic
    odor4pct.executeTrackFunction('fixHTOrientation');
    toc
end

%we are not actually going to stitch tracks any more; I think this causes
%problems by stitching incorrectly past collisions
existsAndDefault('restitchAndTrim',true);
stitchDist = 10; %distance that tracks can be apart to be stitched
stitchFrameDiff = 5; %number of frames that can be between end and next start when stitching
minTrackDist = 75; %minimum distance a track must travel in order to survive trimming
minTrackLength = 400; %after stitching, minimum number of frames a track must have 
minTrackSpeed = 0; %minimum average speed to avoid trimming

if (restitchAndTrim)
    ts = tic;
    %disp('stitching tracks and trimming slow and short ones'); ts = tic(); 
    ntstart = length([odor4pct.expt.track]);
  %  odor4pct.executeExperimentFunction('stitchTracks', stitchFrameDiff, stitchDist);
   % npoststitch = length([odor4pct.expt.track]);
    for j = 1:length(odor4pct.expt)
       meanspeed = odor4pct.expt(j).evaluateTrackExpression('mean(track.getDerivedQuantity(''speed''))');
       dt = sqrt(odor4pct.expt(j).evaluateTrackExpression('max(sum(track.getDerivedQuantity(''displacement'').^2))'));
       valid = [odor4pct.expt(j).track.npts] > minTrackLength & meanspeed > minTrackSpeed & dt > minTrackDist;
       odor4pct.expt(j).track = odor4pct.expt(j).track(valid);
    end
    %analysisRect = [300 125 2275 1875];
    trimRect = [325 150 2250 1850];

    %get rid of tracks that start on the edge
    odor4pct.executeExperimentFunction('pruneTracks', [], trimRect);

    nposttrim = length([odor4pct.expt.track]);
    
    disp(['num tracks start: ' num2str(ntstart) ', after trim: ' num2str(nposttrim)]);
    restitchAndTrim = false;
    toc(ts)
end
existsAndDefault('fixBadTracks', true);
if (fixBadTracks)
    mra = MaggotReAnalyzer();
    mra.contourScale = 3;
    for j = 1:length(odor4pct.expt)
        mraArr{j} = mra.fixExperiment(odor4pct.expt(j),'reextract',false);
    end
    
    disp ('reextracting bad tracks; this may take a while');
    tic
    for k = 1:length(mraArr)
        mra = mraArr{k};
        for j = 1:length(mra)
            mra(j).reExtractTrack([],'resegment',false);
            [mra(j).track.pt.imData] = deal([]);
        end
        disp ([num2str(k) '/' num2str(length(mraArr)) ' completed']);
    end
end
%set segmentation options here
so = MaggotSegmentOptions();
so.curv_cut = .6;

existsAndDefault('resegment', true);
if (resegment)
    %synch up segment options
    [odor4pct.expt.so] = deal(so);
    for j = 1:length(odor4pct.expt)
        [odor4pct.expt(j).track.so] = deal(odor4pct.expt(j).so);
    end
%    odor4pct.executeExperimentFunction('[expt.track.so] = deal(expt.so)');
    
    %recommend running setSegmentSpeeds to let track set its own segment
    %speed based on high curvature regions
    odor4pct.executeTrackFunction('setSegmentSpeeds');
    
    disp('segmenting tracks, this can take a while'); tic
    odor4pct.executeTrackFunction('segmentTrack');
    toc
    resegment = false;
end

varsToClean = {'minLoadPts', 'stitchDist','stitchFrameDiff','minTrackDist = 75','minTrackLength = 400','minTrackSpeed',...
    'ntstart', 'npoststitch', 'meanspeed', 'dt', 'valid', 'nposttrim', 'so', 'varsToClean', 'ts'};
clear(varsToClean{:});

existsAndDefault('saveMe', false);
if (saveMe)
    disp('saving to disk, this can take an ass long time');
    try 
        tic
        save ([basedir 'importedToMatlab.mat'],'odor4pct');
        toc
    catch me
        disp('ok, it choked on that, so I will save individual experiments instead; wait some more');
        tic
        odor4pct.toMatFiles([basedir 'odor4pct']);
        toc
    end
    saveMe = false;
end

fignum = 1;
odor4pct.defaultTitle = 'Ethyl Acetate ~ 80 ppm at midpoint';


figure(fignum); clf(fignum); fignum = fignum+1; 
odor4pct.makeHistogram('theta', deg2rad(-175.5:5:175.5), 'runs', 'r2d',true);
title ('distribution of instantaneous direction');

figure(fignum); clf(fignum); fignum = fignum+1; 

odor4pct.gatherSubField('run', 'startTheta');

