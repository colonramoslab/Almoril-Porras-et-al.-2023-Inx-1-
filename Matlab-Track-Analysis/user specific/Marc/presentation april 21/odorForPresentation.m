savedir = 'C:\Documents and Settings\Marc\My Documents\presentations\';

if (~exist('odordilute', 'var'))
    odordilute = ExperimentSet.fromMatFiles('D:\Marc Processed\maggots\ethyl acetate 4 pct 20 2000\odor4pct');
end
existsAndDefault('trim', true);
if (trim)
    trimrect = [425 130 2150 1850];
    odordilute.executeExperimentFunction('trimTracks', [],trimrect);
end
%%
fignum = 0; 
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
%use to make sure points are selected from a square, not rectangular area
indexpr = 'track.dq.sloc(1,:) > 413 & track.dq.sloc(1,:) < 2163';

binsize = 30;
tx = deg2rad(0:30:330);
[txplot,I] = sort(mod(rad2deg(tx) - 270, 360));
ho = odordilute.makeHistogram('theta', tx, 'run','indsExpression', indexpr,'polar',true);
set(bar (txplot, ho(I)/sum(ho), 'w','LineWidth', 3),'EdgeColor', 'w', 'FaceColor','none','LineWidth',3);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]-15);

xlabel ('heading (degrees)'); ylabel ('fraction of run time'); 
title ('distribution of instantaneous run direction');
embiggen();
saveas(gcf, [savedir 'odor theta distribution.tiff'], 'tiff');
%%
fignum = 2; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
u = get(gca, 'Units');
set(gca, 'Units', 'Pixels');
p = get(gca,'Position');
h = p(4);
w = p(3);
hn = p(3) * 16/20;
dh = h - hn;
p(2) = p(2) + dh/2;
p(4) = hn;
set(gca, 'Position', p)
set(gca, 'Units', u);
t = [eset.expt.track];
t.plotPath('displacement', 'w'); xlim([-800 800]); ylim([-1000 1000]); 
hold on;
plot (0,0,'m.', 'MarkerSize', 50);
set(gca, 'XTick', [], 'YTick', [],'LineWidth', 3)
view(-90,90);
saveas(gcf, [savedir 'odor displacement.tiff'], 'tiff');
%%
fignum = 3; figure(fignum); clf(fignum)
tx = deg2rad(0:45:315);
[txplot,I] = sort(mod(rad2deg(tx) - 270, 360));

set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
h = odordilute.makeReorientationHistogram('theta',tx, 'polar', 'true','minHS',1);
%plot (txplot, h(I), 'w-', 'LineWidth', 3);
set(bar (txplot, h(I), 'LineWidth', 3),'EdgeColor', 'w', 'FaceColor','none','LineWidth',3);
yl = get(gca, 'YLim');
yl(1) = 0;
xlabel ('direction');
ylabel ('reorientation rate (min^{-1})');
set(gca,'XLim', [0 360]-22.5, 'XTick', 0:45:315, 'YLim', yl);
embiggen()
saveas(gcf, [savedir 'odor reorientation rate.tiff'], 'tiff')
%%
tx = deg2rad(0:45:315);
fignum = 4; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
nhs = odordilute.gatherSubField('reorientation', 'numHS');
pd = odordilute.gatherSubField('reorientation', 'prevDir');
nd = odordilute.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([pd;nd]));

[x,meany] = meanyvsx(mod(pd(nhs>0)-pi/2,2*pi), dt(nhs>0), tx);
[x,I] = sort(mod(x - 3*pi/2, 2*pi));
plot(rad2deg(x), rad2deg(meany(I)), 'w-', 'LineWidth', 3);
set(gca,'XLim', [0 360], 'XTick', 0:60:360) 


%%
tx = deg2rad(0:45:315);
fignum = 5; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
nhs = odordilute.gatherSubField('reorientation', 'numHS');
pd = odordilute.gatherSubField('reorientation', 'prevDir');
nd = odordilute.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([pd;nd]));

[x,meany] = meanyvsx(mod(pd(nhs>0)-pi/2,2*pi), abs(dt(nhs>0)), tx);
[x,I] = sort(mod(x - 3*pi/2, 2*pi));
plot(rad2deg(x), rad2deg(meany(I)), 'w-', 'LineWidth', 3);
set(gca,'XLim', [0 360], 'XTick', 0:45:360) 
%%
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

[x,meany] = meanyvsx(mod(pd(nhs>0)-pi/2,2*pi), abs(dt(nhs>0)), tx);

plot(rad2deg(x), rad2deg(meany), 'w-', 'LineWidth', 3);
set(gca,'XLim', [0 360], 'XTick', 0:60:360) 
