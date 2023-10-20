basedir = '\\labnas1\Share\Phototaxis\Data Extracted\Checkerboard\main experiment\';
matbase= '\\labnas1\Share\Phototaxis\Data Extracted\Checkerboard\polys\10000 vs 100\polys_07';
day = [22 24 15 15 23 23 24]; %
time = [1545 1450 1035 1315 1400 1610 1040]; %
datadir = '\\labnas1\Share\Phototaxis\Data\2009\Jul\';

clear fnlist matfn backimfn;
for j = 1:length(day)
    fnlist{j} = [basedir num2str(day(j)) '_' num2str(time(j)) '_tracks.bin'];
    matfn{j} = [matbase num2str(day(j)) '09_' num2str(time(j)) '.mat'];
    backimfn{j} = [datadir num2str(day(j)) '\' num2str(time(j)) '\' num2str(time(j)) '_background\' num2str(time(j)) '_background_0.jpg']; 
end


if (~exist('eset', 'var'))
    eset = ExperimentSet.fromFiles(fnlist{:},'minpts',100);
end
if (~exist('checker', 'var'))
    checker = ExperimentSet();
    for j = 1:length(eset.expt)
        expt2(j) = CheckerExperiment(eset.expt(j), matfn{j});
    end
    checker.expt = expt2;
    checker.executeExperimentFunction('addSpatialInfo');
end

ecl = ESetCleaner;
ecl.minSpeed = 1.8;
ecl.minPts = 500;
ecl.minHTValid = 0.98; 
existsAndDefault ('clean', true);
if (clean)
    ecl.clean(checker);
    clean = false;
end
existsAndDefault ('fixht', true);
if (fixht)
    checker.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

existsAndDefault ('segment', true);
if (segment)
    checker.executeTrackFunction('setSegmentSpeeds');
    checker.executeTrackFunction('segmentTrack');
    segment = false;
end



return
%%
close all; fignum = 0;

fignum = fignum + 1; figure(fignum); clf(fignum)
h1 = checker.makeReorientationHistogram('boundaryDistance', -100:5:100);
h2 = checker.makeReorientationHistogram('headBoundaryDistance', -100:5:100);
plot (-100:5:100, h1, -100:5:100, h2); title ('reorienation rate vs. distance from boundary'); legend ('mid', 'head');

fignum = fignum + 1; figure(fignum); clf(fignum)
make2DReorientationHistogram(checker, 'headBoundaryDistance', -45:3:45, 'boundaryDistance', -45:3:45);
axis square; axis equal; shading interp
hold on
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
plot (xl, [0 0], 'k--',[0 0],yl, 'k--', 'LineWidth', 2);
plot (-45:3:45, -45:3:45, 'm--','LineWidth', 2);



%%
mbd = checker.gatherFromSubField('reorientation', 'boundaryDistance', 'position', 'start');
hbd = checker.gatherFromSubField('reorientation', 'headBoundaryDistance', 'position', 'start');
mbdall = checker.gatherFromSubField('run', 'boundaryDistance', 'indsExpression', 'notlast');
hbdall = checker.gatherFromSubField('run', 'headBoundaryDistance', 'indsExpression', 'notlast');

%%
xaxis = [-100:5:100];
yaxis = [-20:2:20];
imreo = makeIm(mbd, hbd-mbd, xaxis, yaxis);
imrun = makeIm(mbdall, hbdall-mbdall, xaxis, yaxis);


figure(11);
pcolor (xaxis,yaxis,240*imreo./imrun); shading interp

figure(12);

plot (mbdall, hbdall-mbdall, 'y.', mbd, hbd-mbd, 'b*',[-100 100],[-100 100],'r-', [-100 100],[0 0], 'g-', [0 0],[-100 100], 'g-'); axis ([-20 20 -20 20])

xlabel ('mid position (>0 in dark)'); ylabel ('head position');
%%
figure(1);
tx = deg2rad(-180:45:180);
h1 = checker.makeReorientationHistogram ('thetaToBound', tx, 'r2d', true, 'minHS', 1, 'validname', 'onboundary');
h2 = checker.makeReorientationHistogram ('thetaToBound', tx, 'r2d', true, 'validname', 'onboundary');
h3 = checker.makeReorientationHistogram ('thetaToBound', tx, 'r2d', true, 'maxHS', 0, 'validname', 'onboundary');
plot (rad2deg(tx), h1, rad2deg(tx), h2, rad2deg(tx), h3); title ('reo-rate'); ylabel ('min^{-1}'); xlabel ('degrees'); legend ('> 0 headsweeps', 'all reos', '0 headsweeps');

figure(2);
plot (rad2deg(tx), h3./h2); title ('fraction of reorientations that are pauses');
%%
figure(3);
rd = checker.gatherFromSubField('reorientation', 'thetaToBound', 'position', 'start');
onb = logical(checker.gatherFromSubField('reorientation', 'onboundary', 'position', 'start'));
nhs = checker.gatherSubField('reorientation', 'numHS');

[x,meany] = meanyvsx(rd(onb), nhs(onb), tx); 
plot(rad2deg(x), meany); title ('num HS/reorientation vs. start direction'); xlabel ('start direction (degrees)'); ylabel ('mean HS/reo');

%%
figure(4);
tx = deg2rad(-180:30:180);
%sd = checker.gatherFromSubField('reorientation', 'theta', 'position', 'start');
sd = checker.gatherSubField('reorientation', 'prevDir');
ed = checker.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([sd;ed]));
[x,meany] = meanyvsx(rd(onb), rad2deg(dt(onb)), tx); 
plot(rad2deg(x), meany); title ('mean angle change vs. angle to boundary'); xlabel ('start direction (degrees)'); ylabel ('mean HS/reo');

figure(5)
[x,meany] = meanyvsx(rd(onb), rad2deg(abs(dt(onb))), tx); 
plot(rad2deg(x), meany); title ('mean magnitude angle change vs. angle to boundary'); xlabel ('start direction (degrees)'); ylabel ('mean HS/reo');

%%
tx = deg2rad(-120:20:120);
hstd = checker.gatherSubField('headSwing', 'tailDir');
hshd = checker.gatherSubField('headSwing', 'headDir');
hsacc = checker.gatherSubField('headSwing', 'accepted');
bounddir = checker.gatherFromSubField('headSwing', 'boundarytheta', 'position', 'mean');
onb = logical(checker.gatherFromSubField('headSwing', 'onboundary', 'position', 'start'));


hstd = diff(unwrap([bounddir;hstd]));
hshd = diff(unwrap([bounddir;hshd]));

hschange = abs(hshd)-abs(hstd);
hschange(sign(hshd) ~= sign(hstd)) = NaN;
hschangeall = hschange;
figure(6)
[x,meany] = meanyvsx(hschange(onb), hsacc, tx); 
plot (rad2deg(x), meany); title ('probability of acceptance vs. change w.r.t boundary')
sum(hschangeall < 0 & hsacc) / sum(hschangeall < 0)
sum(hschangeall > 0 & hsacc) / sum(hschangeall > 0)

%%
tx = deg2rad(-120:20:120);
hstd = checker.gatherSubField('firsths', 'tailDir');
hshd = checker.gatherSubField('firsths', 'headDir');
hsacc = checker.gatherSubField('firsths', 'accepted');
bounddir = checker.gatherFromSubField('firsths', 'boundarytheta', 'position', 'mean');
onb = logical(checker.gatherFromSubField('firsths', 'onboundary', 'position', 'start'));


hstd = diff(unwrap([bounddir;hstd]));
hshd = diff(unwrap([bounddir;hshd]));

hschange = abs(hshd)-abs(hstd);
hschange(sign(hshd) ~= sign(hstd)) = NaN;
figure(8)
plot(rad2deg(tx), hist(hschange,tx), rad2deg(tx), hist(hschangeall,tx), rad2deg(tx),hist(hschangeall,tx)- hist(hschange,tx));
title ('histogram of head sweeps that do not cross 0 or 180'); xlabel ('change wrt. boundary'); ylabel ('num HS');
legend ('first hs', 'all hs', 'all - first hs');
%%
figure(7)
alltowarddark = sum(hschangeall < 0 & abs(hschangeall) < 2*pi/3)/sum(abs(hschangeall < 2*pi/3));
firsttowarddark =  sum(hschange < 0 & abs(hschange) < 2*pi/3)/sum(abs(hschange < 2*pi/3));
otherstowarddark = (sum(hschangeall < 0 & abs(hschangeall) < 2*pi/3) - sum(hschange < 0 & abs(hschange) < 2*pi/3)) / ...
    (sum(abs(hschangeall < 2*pi/3)) - sum(abs(hschange < 2*pi/3)));

alltowarddark
firsttowarddark
otherstowarddark

alltowardlight = sum(hschangeall > 0 & abs(hschangeall) < 2*pi/3)/sum(abs(hschangeall < 2*pi/3));
firsttowardlight =  sum(hschange > 0 & abs(hschange) < 2*pi/3)/sum(abs(hschange < 2*pi/3));
otherstowardlight = (sum(hschangeall > 0 & abs(hschangeall) < 2*pi/3) - sum(hschange > 0 & abs(hschange) < 2*pi/3)) / ...
    (sum(abs(hschangeall < 2*pi/3)) - sum(abs(hschange < 2*pi/3)));

alltowardlight
firsttowardlight
otherstowardlight
bar ([alltowarddark, firsttowarddark,otherstowarddark]);
%set(gca, 'XTickLabel', {'all headswings', 'first headswing', 'other headswings'})


