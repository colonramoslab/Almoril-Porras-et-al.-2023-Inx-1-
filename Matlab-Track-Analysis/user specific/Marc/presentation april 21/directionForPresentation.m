basedir = '\\labnas1\Share\Phototaxis\Data Extracted\Direction\45degrees\';
savedir = 'C:\Documents and Settings\Marc\My Documents\presentations\';

mmPerPixel = 0.096;
if (~exist('eset', 'var'))
    eset = ExperimentSet.fromFiles(basedir,'minpts',100);
end

existsAndDefault('dotrim',true);
if (dotrim)
    trimRect = [360 75 2260 1875];
    eset.executeExperimentFunction('trimTracks', [], trimRect);
    dotrim = false;
end

ecl = ESetCleaner;
ecl.minSpeed = 2;
ecl.minPts = 500;
ecl.minHTValid = 0.98; 
ecl.minDist = 50;
ecl.askFirst = false;
existsAndDefault ('clean', true);
if (clean)
     ecl.clean(eset);
    clean = false;
end
existsAndDefault ('fixht', true);
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

existsAndDefault ('dosegment', true);
if (dosegment)
    eset.executeTrackFunction('setSegmentSpeeds');
    eset.executeTrackFunction('segmentTrack');
    dosegment = false;
end
%%
fignum = 0; 
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

t = [eset.expt.track];
t.plotPath('sloc', 'w'); xlim([500 2000]); axis equal;  axis equal 
view(90,90);
x0 = 500;
y0 = 75;
xtl = 0:1:15;
ytl = 0:1:16;
p = get(gca,'Position');
h = p(4);
w = p(3);
hn = p(3) * 15/18;
dh = h - hn;
p(2) = p(2) + dh/2;
p(4) = hn;
set(gca, 'Position', p)
set(gca,'XTick',10*xtl/mmPerPixel + x0, 'XTickLabel', xtl, 'YTick',10*ytl/mmPerPixel + y0, 'YTickLabel', ytl); 
xlim([500 2000]); ylim([75 1875]);
xlabel('cm'); ylabel('cm')
embiggen();
%saveas(gcf, [savedir 'direction tracks.tiff'], 'tiff')
%%
fignum = 2; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

t([t.startFrame] < 1000).plotPath('displacement', 'w'); xlim([-750 750]); ylim([-900 900]); hold on
plot (0,0,'c.','MarkerSize',40)
view(90,90);
%x0 = 500;
%y0 = 75;
xtl = -8:1:8;
ytl = -10:1:10;
u = get(gca, 'Units');
set(gca, 'Units', 'Pixels'); 
p = get(gca,'Position');
h = p(4);
w = p(3);
hn = p(3) * 15/18;
dh = h - hn;
p(2) = p(2) + dh/2;
p(4) = hn;
set(gca, 'Position', p)
set(gca, 'Units', u);
set(gca,'XTick',10*xtl/mmPerPixel, 'XTickLabel', xtl, 'YTick',10*ytl/mmPerPixel, 'YTickLabel', ytl); 
%xlim([500 2000]); ylim([75 1875]);
xlabel('cm'); ylabel('cm')
embiggen();
saveas(gcf, [savedir 'direction displacement.tiff'], 'tiff')

%%
tx = deg2rad(0:30:330);
fignum = 3; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
h = eset.makeHistogram('theta',tx, 'polar', 'true');
[tx2,I] = sort(mod(rad2deg(tx)-90,360));

set(bar (tx2, h(I)/sum(h), 'w','LineWidth', 3),'EdgeColor', 'w', 'FaceColor','none','LineWidth',3);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]-15);

xlabel ('heading (degrees)'); ylabel ('fraction of run time'); 
title ('distribution of instantaneous direction');
%{
plot (tx2, h(I)/sum(h), 'w-', 'LineWidth', 3);
yl = get(gca, 'YLim');
yl(1) = 0;
xlabel ('direction');
ylabel ('fraction of run time');
set(gca,'XLim', [0 360], 'XTick', 0:60:360, 'YLim', yl);
%}
embiggen()
saveas(gcf, [savedir 'direction run direction.tiff'], 'tiff')

%%
tx = deg2rad(0:30:330);
fignum = 4; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
h = eset.makeReorientationHistogram('theta',tx, 'polar', 'true');
[tx2,I] = sort(mod(rad2deg(tx)-90,360));

set(bar (tx2, h(I)/sum(h), 'w','LineWidth', 3),'EdgeColor', 'w', 'FaceColor','none','LineWidth',3);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]-15);

xlabel ('heading (degrees)'); ylabel ('reorientation rate'); 
title ('reorientation rate vs. direction');
%{
plot (tx2, h(I), 'w-', 'LineWidth', 3);
yl = get(gca, 'YLim');
yl(1) = 0;
xlabel ('direction');
ylabel ('reorientation rate (min^{-1})');
set(gca,'XLim', [0 360], 'XTick', 0:60:360, 'YLim', yl);
%}
embiggen()
saveas(gcf, [savedir 'direction reorientation rate.tiff'], 'tiff')
%%
tx = deg2rad(0:30:330);
fignum = 5; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
hsd = mod(eset.gatherSubField('firsths', 'headDir') + pi/12, 2*pi)-pi/12;
acc = logical(eset.gatherSubField('firsths', 'accepted'));
[tx2,I] = sort(mod(rad2deg(tx)-90,360));
h = hist(hsd,tx);
h2 = hist(hsd(acc), tx);
set(bar (tx2, h(I)/sum(h), 'w','LineWidth', 3),'EdgeColor', 'w', 'FaceColor','none','LineWidth',3);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'XTick', 0:60:360,'XLim',[0 360]-15);
title ('first head sweep direction');
ylabel ('fraction of head sweeps');
xlabel ('direction');
embiggen();
saveas(gcf, [savedir 'direction head sweep bias.tiff'], 'tiff')

%%
dtx = -180:10:180;
fignum = 5; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
rs = eset.gatherSubField('run','startTheta');
re = eset.gatherSubField('run','endTheta');
dt = rad2deg(diff(unwrap([rs;re])));
to0 = find(abs(rs - pi/2) < pi/4);
to180 = find(abs(rs + pi/2) < pi/4);
ton90 = find(abs(rs) < pi/4);
top90 = find(abs(rs) > 3*pi/4);
subplot(2,2,1);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
bar(dtx, hist(dt(to0), dtx), 'r');
title ({'run starts towards 0',['mean = ' num2str(mean(dt(to0))) '; stdev = ' num2str(std(dt(to0)))]});
embiggen();
subplot(2,2,2);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

bar(dtx, hist(dt(to180), dtx), 'c');
title ({'run starts towards 180',['mean = ' num2str(mean(dt(to180))) '; stdev = ' num2str(std(dt(to180)))]});
embiggen();
subplot(2,2,3);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

bar(dtx, hist(dt(ton90), dtx), 'g');
title ({'run starts towards -90',['mean = ' num2str(mean(dt(ton90))) '; stdev = ' num2str(std(dt(ton90)))]});
embiggen();
subplot(2,2,4);
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

bar(dtx, hist(dt(top90), dtx), 'y');
title ({'run starts towards +90',['mean = ' num2str(mean(dt(top90))) '; stdev = ' num2str(std(dt(top90)))]});
embiggen();
saveas(gcf, [savedir 'direction run steering.tiff'], 'tiff')
%%
tx = deg2rad(0:30:330);
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');
nhs = eset.gatherSubField('reorientation', 'numHS');
pd = eset.gatherSubField('reorientation', 'prevDir');
nd = eset.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([pd;nd]));

[x,meany] = meanyvsx(mod(pd(nhs>0)-pi/2,2*pi), dt(nhs>0), tx);

plot(rad2deg(x), rad2deg(meany), 'w-', 'LineWidth', 3);
set(gca,'XLim', [0 360], 'XTick', 0:60:360) 
%%
fignum = fignum + 1; figure(fignum); clf(fignum)
set(gcf, 'Color', 'k','InvertHardCopy', 'off');
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1],'LineWidth',2,...
    'DefaulttextColor', [1 1 1], 'NextPlot', 'replaceChildren','box','on');

[x,meany] = meanyvsx(mod(pd(nhs>0)-pi/2,2*pi), abs(dt(nhs>0)), tx);

plot(rad2deg(x), rad2deg(meany), 'w-', 'LineWidth', 3);
set(gca,'XLim', [0 360], 'XTick', 0:60:360) 
