basedir = '\\labnas1\Share\Phototaxis\Data Extracted\Checkerboard\main experiment\';
matbase= '\\labnas1\Share\Phototaxis\Data Extracted\Checkerboard\polys\10000 vs 100\polys_07';
day = [22 24 15 15 23 23 24]; %
time = [1545 1450 1035 1315 1400 1610 1040]; %
datadir = '\\labnas1\Share\Phototaxis\Data\2009\Jul\';
savedir = 'C:\Documents and Settings\Marc\My Documents\presentations\';
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

existsAndDefault ('dosegment', true);
if (dosegment)
    checker.executeTrackFunction('setSegmentSpeeds');
    checker.executeTrackFunction('segmentTrack');
    dosegment = false;
end





%%
existsAndDefault('saveFigs', false);
close all; fignum = 0;
mmperpixel = 0.066; %CAN'T BE RIGHT
mmperpixel = 0.086; %Maybe she said 86 microns / pixel?
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
dx = -130:5:130;
h1 = checker.makeReorientationHistogram('boundaryDistance', dx);
h2 = checker.makeReorientationHistogram('headBoundaryDistance', dx);
plot (dx*mmperpixel, h1,'w', 'LineWidth', 3); 
set(title ('reorienation rate vs. distance from boundary'), 'Color', 'w'); 
%set(legend ('head', 'midpoint'),'TextColor','w');
xlabel ('dist (mm)');
ylabel ('rate (min^{-1})');
embiggen();

fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
plot (dx*mmperpixel, h1,'w', 'LineWidth', 3); 
hold on;
yl = get(gca,'YLim');
patch([-30 30 30 -30]*mmperpixel, [yl(1) yl(1) yl(2) yl(2)], 'c', 'FaceAlpha', 0.3);
set(title ('reorienation rate vs. distance from boundary'), 'Color', 'w'); 
%set(legend ('head', 'midpoint'),'TextColor','w');
xlabel ('dist (mm)');
ylabel ('rate (min^{-1})');
embiggen();

dx = -45:3:45;
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
h = make2DReorientationHistogram(checker, 'boundaryDistance', dx, 'headBoundaryDistance', dx);
pcolor (dx*mmperpixel,dx*mmperpixel,h);
axis square; axis equal; shading interp;
colormap hot
colorbar vert
hold on
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
%plot (xl, [0 0], 'w--',[0 0],yl, 'w--', 'LineWidth', 2);
plot (-45:3:45, -45:3:45, 'w--','LineWidth', 2);
xlabel ('mid point position');ylabel('head position');
%%
fignum = 6; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
tx = deg2rad(0:30:330);
h2 = checker.makeReorientationHistogram ('thetaToBound', tx, 'r2d', true, 'validname', 'onboundary','polar',true);
plot (rad2deg(tx), h2, 'w-','LineWidth',3); 
set(title ('Reorientation Rate vs. Angle of Velocity To Boundary'),'Color','w');
ylabel ('rate min^{-1}'); xlabel ('heading');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:330,'XLim',[0 330]);
embiggen()

if(saveFigs)
    saveas(gcf, [savedir 'reorientation rate checker linear.tiff'], 'tiff');
end

fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
ax = newplot();
set(ax, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','FontSize',14);
set(polar (ax, [tx tx(1)], [h2 h2(1)],'w'),'LineWidth',3)
set(title('reorientation rate vs. angle to boundary'),'Color','w','FontSize',16);
return
%%

tx = deg2rad(0:30:360);
rd = mod(checker.gatherFromSubField('reorientation', 'thetaToBound', 'position', 'start'),2*pi);
onb = logical(checker.gatherFromSubField('reorientation', 'onboundary', 'position', 'start'));
nhs = checker.gatherSubField('reorientation', 'numHS');
%{
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

[x,meany] = meanyvsx(rd(onb), nhs(onb), tx); 
plot(rad2deg(x), meany,'w'); 
set(title ('num HS/reorientation vs. start direction'),'TextColor','w');
xlabel ('start direction (degrees)'); 
ylabel ('mean HS/reo');
%}

fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
%tx = deg2rad(-180:30:180);
%sd = checker.gatherFromSubField('reorientation', 'theta', 'position', 'start');
sd = checker.gatherSubField('reorientation', 'prevDir');
ed = checker.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([sd;ed]));
[x,meany] = meanyvsx(mod(rd(onb),2*pi), rad2deg(dt(onb)), tx); 
reox = x;
reomeany=meany;
plot(rad2deg(x), meany,'w','LineWidth',3); 
set(title ('mean angle change vs. angle to boundary'),'Color','w'); 
xlabel ('start direction (degrees)'); 
ylabel ('mean angle change');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]);

embiggen();
if saveFigs
    saveas(gcf, [savedir 'reorientation direction at boundary checker.tiff'],'tiff');
end
%%
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
[x,meany] = meanyvsx(rd(onb), rad2deg(abs(dt(onb))), tx); 
plot(rad2deg(x),  meany,'w','LineWidth',3); 
set(title ('magnitude of angle change vs. angle to boundary'),'Color','w'); 
xlabel ('start direction (degrees)'); 
ylabel ('magnitude of angle change');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]);

embiggen();
if saveFigs
    saveas(gcf, [savedir 'reorientation magnitude at boundary checker.tiff'],'tiff');
end
%%
tx = deg2rad(0:30:360);

fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
rd = mod(checker.gatherFromSubField('firsths', 'thetaToBound', 'position', 'start'),2*pi);
onb = logical(checker.gatherFromSubField('firsths', 'onboundary', 'position', 'start'));

dt = checker.gatherSubField('firsths', 'maxTheta');
[x,meany] = meanyvsx(mod(rd(onb),2*pi), rad2deg(dt(onb)), tx); 
plot(rad2deg(reox), reomeany, 'c--', 'LineWidth',2); hold on
plot(rad2deg(x), meany,'w','LineWidth',3); 
set(title ('mean first headsweep angle change vs. angle to boundary'),'Color','w'); 
xlabel ('start direction (degrees)'); 
ylabel ('mean angle change');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]);
embiggen();
set(legend ('reorientation angle', 'first headsweep angle','Location','NorthWest'),'TextColor','w','FontSize',14);
if saveFigs
    saveas(gcf, [savedir 'headsweep direction at boundary checker.tiff'],'tiff');
end
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


