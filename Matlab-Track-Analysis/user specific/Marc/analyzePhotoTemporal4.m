basedir ='\\labnas1\Share\Ashley Extracted\Temporal\redo Temporal2\old 60s';
if (~exist('photo60', 'var'))
    %minpts = 500;
    photo60 = ExperimentSet.fromFiles(basedir,'minpts',50);
end
if (~isa(photo60.expt(1), 'TemporalLightExperiment'))
    expt = photo60.expt;
    for j = 1:length(photo60.expt)
        expt2(j) = TemporalLightExperiment(expt(j), 240);
    end
    photo60.expt = expt2;
end

%%
existsAndDefault('stitchTracks', true);
if (stitchTracks)
    %first auto-stitch without approval any frame diff of 1 with distance
    %less than 10 
    %photo60.executeExperimentFunction('stitchTracks', 1, 10, 'interactive', false);
    
    %next autostitch without approval any frame diff of 2 with distance < 4
    photo60.executeExperimentFunction('stitchTracks', 2, 4, 'interactive', false);
    
    frameDiff = 4;
    maxDist = 20;
    photo60.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', true);
    stitchTracks = false;
end
close all;
ecl = ESetCleaner();
ecl.minHTValid = 0;
ecl.minPts = 500;
ecl.askFirst = true;
ecl.clean(photo60);
close all;
photo60.defaultTitle = 'Temporal Phototaxis 60 sec cycle';
ecl = ESetCleaner();
ecl.minSpeed = 2.7;
ecl.minDist = 100;
ecl.minPts = 500;
ecl.clean(photo60);
close all;
%photo60.makeHistogram('speed', 0:0.25:7);
%%
fignum = 0;
existsAndDefault('makeSpeedFigs', true);

if (makeSpeedFigs)
    tx = 0:2:240;
    fignum = fignum + 1; figure(fignum); clf(fignum);
    photo60.meanField2vsField1('timeon', 'speed', tx); hold on;
    %yl = get(gca, 'Ylim'); yl(1) = 0;
    yl = [0 5];
    plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);

    fignum = fignum + 1; figure(fignum); clf(fignum);
    photo60.meanField2vsField1('timeoff', 'speed', tx); hold on;
    %yl = get(gca, 'Ylim'); yl(1) = 0;
    yl = [0 5];
    plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);

   
end
%%

existsAndDefault('fixHT', true);
if (fixHT)
    disp('fixing HT orientation');
    photo60.executeTrackFunction('fixHTOrientation','mintime',10);
    fixHT = false;
end
existsAndDefault('dosegment', 'true');
if (dosegment)
    disp('setting segment params automatically');
    photo60.executeTrackFunction('setSegmentSpeeds');
    disp('segmenting');
    photo60.executeTrackFunction('segmentTrack');
    dosegment = false;
end
%%
fignum = fignum + 1; figure(fignum); clf(fignum);
r1 = photo60.makeReorientationHistogram('timeon', tx);
r2 = photo60.makeReorientationHistogram('timeon', tx, 'minHS', 1);
plot (tx, r1, tx, r2); legend ('all', 'no pauses');
title ('reorientation rate vs. time on');
xlabel('time since light turned on (s)'); ylabel ('rate');
fignum = fignum + 1; figure(fignum); clf(fignum);
photo60.makeReorientationHistogram('timeoff', tx);

%%
fignum = fignum + 1; figure(fignum); clf(fignum);
binsize = 10;
tx2 = (binsize/2):binsize:(240-binsize/2);
hst = photo60.makeHistogram('timeon', tx, 'hsstart');
hst2 = photo60.makeHistogram('timeon', tx, 'reo','start');
allt = photo60.makeHistogram('timeon', tx, 'run', 'notlast');
plot (tx, hst./allt * 240, tx, hst2./allt * 240); xlabel ('light time'); ylabel ('rate (min^{-1})'); title ('headsweep initiation rate vs. cycle');
embiggen();

fignum = fignum+1; figure(fignum); clf(fignum);
photo60.meanField2vsField1('timeoff', 'speed', tx, 'runs'); hold on;
plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);
title ('speed in runs only'); hold off

%%
tx2 = 0:5:240;
fignum = fignum+1; figure(fignum); clf(fignum);
nhs = photo60.gatherSubField('reorientation', 'numHS');
reot = photo60.gatherFromSubField('reorientation', 'timeon', 'position', 'start');
[x,meany] = meanyvsx (reot, nhs, tx2);
[x2,meany2] = meanyvsx(reot(nhs > 0), nhs(nhs > 0), tx2);
plot (x,meany, x2, meany2); title ('headswings / reorientation'); legend ('all reos', 'omit pauses');
xlabel ('time on');
ylabel ('hs/reo')