if (~exist ('checkeset', 'var'))
    checkeset = ExperimentSet.fromMatFiles('E:\phototaxis\checkerboard\matfiles\cs_checker_3cm');
    checkeset.executeTrackFunction('setSegmentSpeeds');
    checkeset.executeTrackFunction('segmentTrack');
end
if (~exist('cec', 'var'))
     load E:\phototaxis\checkerboard\checkcalc.mat
end
figure(1); clf(1);
pcolor (cec.rx, cec.ry, cec.morphedCheckerIm); shading flat; colormap gray; ylim([1 20]); axis equal; ylim([1 20]); hold on
checkeset.executeTrackFunction('plotPath','sloc', 'm-'); set(gca, 'XTick', [], 'YTick', []);

if (~exist ('ccalc', 'var'))
    ccalc = checkerCalculations (checkeset);
end

if (~exist ('adcheck', 'var'))
    cno = checkerNavigationAnalysis;
    adcheck = checkerNavigationAnalysis(checkeset, cno, ccalc);
end

hset = checkerNavigationFigures(adcheck);