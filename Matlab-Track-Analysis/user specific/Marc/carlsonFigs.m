%generate figures for presentation to john carlson 
%compare odor guided navigation to control
%March 06, 2010
savedir = 'D:\Marc Processed\CarlsonFigsMar2010\';


if (~exist ('control', 'var'))
    disp ('loading control files analyzed'); 
    tic
    load 'D:\Marc Processed\maggots\control\importedToMatlab.mat'
    toc
end
if (~exist ('odordilute', 'var'))
    disp ('loading control files analyzed');
    tic
    load 'D:\Marc Processed\maggots\ethyl acetate 4 pct 20 2000\importedToMatlab.mat'
    toc
end

%use to make sure points are selected from a square, not rectangular area
indexpr = 'track.dq.sloc(1,:) > 413 & track.dq.sloc(1,:) < 2163';

fignum = 0;
fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
binsize = 30;
tx = thetaAxis(binsize, 'd', 'r');
txd = thetaAxis(binsize, 'd', 'd');
hc = control.makeHistogram('theta', tx, 'run','indsExpression', indexpr);
sum(hc);
ho = odordilute.makeHistogram('theta', tx, 'run','indsExpression', indexpr);
plot (txd, hc/sum(hc), 'b-', txd, ho/sum(ho), 'r-', 'LineWidth', 3);
xlabel ('heading (degrees)'); ylabel ('fraction of run time'); title ('distribution of instantaneous run direction');
legend ('control', 'odor gradient')
emsmallen(gca,'FontSize',8);

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
txd = [-180:90:180];
tx = deg2rad(txd);
hc = control.makeReorientationHistogram('theta', tx,'indsExpression', 'track.dq.eti > 60');
ho = odordilute.makeReorientationHistogram('theta', tx, 'indsExpression','track.dq.eti > 60');
hc(1) = mean([hc(1) hc(end)]);
hc = hc(1:(end-1));
ho(1) = mean([ho(1) ho(end)]);
ho = ho(1:(end-1));
txd = txd(1:(end-1));
bar (txd, [hc;ho]');
xlabel ('heading (degrees)'); ylabel ('reorientation rate (min^{-1})'); title ('reorientation rate vs. instantaneous run direction');
legend ('control', 'odor gradient')
xtl = get(gca, 'XTickLabel');
xtl(:,5:6) = repmat(' ', 4, 2);
xtl(1,:) = '+/-180';
set(gca, 'XTickLabel', xtl);
emsmallen(gca,'FontSize',8);

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
txd = [-180:90:180];
tx = deg2rad(txd);
hc = control.makeHistogram('theta', tx, 'runstart','indsExpression', 'track.dq.eti > 60');
ho = odordilute.makeHistogram('theta', tx, 'runstart','indsExpression', 'track.dq.eti > 60');
hc(1) = sum([hc(1) hc(end)]);
hc = hc(1:(end-1));
ho(1) = sum([ho(1) ho(end)]);
ho = ho(1:(end-1));
nhc = sum(hc);
nho = sum(ho);
pc = hc/sum(hc);
po = ho/sum(ho);
ec = sqrt(nhc .* pc .*(1 -pc))/nhc;
eo = sqrt(nho .* po .*(1-po))/nhc;

txd = txd(1:(end-1));
errorbar (txd, pc, ec, 'b.', 'MarkerSize', 30, 'LineWidth', 3); hold on
errorbar (txd, po, eo, 'r.', 'MarkerSize', 30, 'LineWidth', 3); hold off
title ({'initial run direction histogram','errorbars represent uncertainty due to counting statistics only'});
xlabel ('intial run heading');
ylabel ('fraction of runs');
legend ('control', 'odor gradient')
set (gca, 'XTick', [-180 -90 0 90]);
xtl = get(gca, 'XTickLabel');
xtl2 = repmat (' ', 4, 6);
xtl2(:,3:6) = xtl;
xtl = xtl2;
xtl(1,:) = '+/-180';
set(gca, 'XTickLabel', xtl);
emsmallen(gca,'FontSize',8);

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
etx = 0:30:900;
hc = control.makeReorientationHistogram('eti', etx);
ho = odordilute.makeReorientationHistogram('eti', etx);
plot (etx, hc, 'b-', etx, ho, 'r-', 'LineWidth', 3); xlim([0 800]);
title ('reorientation rate vs. time');
xlabel('s'); ylabel ('reorientation rate (min^{-1})'); 

legend ('control', 'odor gradient')
emsmallen(gca,'FontSize',8);

%{
fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
etx = 0:60:900;
[ce,cv] = control.meanField2vsField1('eti','vel',etx,'runs');
[oe,ov] = odordilute.meanField2vsField1('eti','vel',etx,'runs');
plot (ce,cv(2,:),'b-',ce,cv(1,:),'b--', oe, ov(2,:), 'r-', oe, ov(1,:), 'r--','LineWidth',2);
title ('drift velocity vs. time');
legend('control y', 'control x', 'odor y', 'odor x'); 
emsmallen(gca,'FontSize',8);
%}
fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
t = [control.expt.track];
%analysisRect = [300 125 2275 1875];
%trimRect = [325 150 2250 1850];

for j = 1:length(t)
    clocend(:,j) = t(j).dq.sloc(:,end);
    ctimeend(j) = t(j).dq.eti(end);
end
t = [odordilute.expt.track];
%analysisRect = [300 125 2275 1875];
%trimRect = [325 150 2250 1850];

for j = 1:length(t)
    olocend(:,j) = t(j).dq.sloc(:,end);
    otimeend(j) = t(j).dq.eti(end);
end

ctop = find(clocend(2,:) > 1850);
cbottom = find(clocend(2,:) < 150);
cright = find(clocend(1,:) > 2250);
cleft = find(clocend(1,:) < 325);
cedge = [ctop, cright, cleft];

otop = find(olocend(2,:) > 1850);
obottom = find(olocend(2,:) < 150);
oright = find(olocend(1,:) > 2250);
oleft = find(olocend(1,:) < 325);
oedge = [otop, oright, oleft];

ctime = ctimeend(cedge);
otime = otimeend(oedge);

ctime = ctime(ctime < 900);
otime = otime(otime < 900);

cended = zeros(1,900);
cended(round(ctime)) = 1;
cended = cumsum(cended);
oended = zeros(1,900);
oended(round(otime)) = 1;
oended = cumsum(oended);

ctime = ctimeend(cbottom);
otime = otimeend(obottom);

ctime = ctime(ctime < 900);
otime = otime(otime < 900);

cendedb = zeros(1,900);
cendedb(round(ctime)) = 1;
cendedb = cumsum(cendedb);
oendedb = zeros(1,900);
oendedb(round(otime)) = 1;
oendedb = cumsum(oendedb);


plot (1:900, cended/length([control.expt.track]), 'b--', 1:900, oended/length([odordilute.expt.track]), 'r--',...
    1:900, cendedb/length([control.expt.track]), 'b-', 1:900, oendedb/length([odordilute.expt.track]), 'r-', 'LineWidth', 3);
xlabel('elapsed time');
ylabel('approximate fraction of tracks that have hit the edge');
title ('animals lost to edge of plate');
legend('control all sides but bottom', 'odor all sides but bottom', 'control bottom side', 'odor bottom side','Location','Best');
emsmallen(gca,'FontSize',8);
existsAndDefault('savefigs', false);

%cd = control.gatherField('displacement');
%ce = control.gatherField('eti');
[ce,cmd] = control.meanField2vsField1('eti', 'displacement', 0:30:900);

%od = odordilute.gatherField('displacement');
%oe = odordilute.gatherField('eti');
[oe,omd] = odordilute.meanField2vsField1('eti', 'displacement', 0:30:900);

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
hold on
odordilute.executeTrackFunction('plotPath', 'sloc', 'k-', 'LineWidth', 0.3);
title ('Raw Paths - Odor');
axis equal
emsmallen(gca,'FontSize',8);

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
hold on
control.executeTrackFunction('plotPath', 'sloc', 'k-', 'LineWidth', 0.3);
axis equal
title ('Raw Paths - Control');
emsmallen(gca,'FontSize',8);


fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
hold on
control.executeTrackFunction('plotPath', 'displacement', 'k-', 'LineWidth', 0.3,'Color', [0.5 0.5 0.5]);
plot (cmd(1,:), cmd(2,:), 'r-', 'LineWidth', 2);
plot (0,0,'k.', 'MarkerSize', 10);
axis([-1250 1250 -1250 1250]);
axis equal
title ('Displaced Paths and Mean Displacement - Control');
emsmallen(gca,'FontSize',8);


fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
hold on
odordilute.executeTrackFunction('plotPath', 'displacement', 'k-', 'LineWidth', 0.3,'Color', [0.5 0.5 0.5]);
plot (omd(1,:), omd(2,:), 'r-', 'LineWidth', 2);
plot (0,0,'k.', 'MarkerSize', 10);
axis([-1250 1250 -1250 1250]);
axis equal
title ('Displaced Paths and Mean Displacement - Odor');
emsmallen(gca,'FontSize',8);


if (savefigs)
    for j = 1:fignum
        saveas(j, [savedir 'figure' num2str(j) '.pdf']);
        saveas(j, [savedir 'figure' num2str(j) '.eps'], 'eps2c');
    end
    savefigs = false;
end



ctd = control.gatherSubField('headSwing', 'tailDir');
cta = control.gatherSubField('headSwing', 'accepted');
ctt = control.gatherSubField('headSwing', 'maxTheta');

otd = odordilute.gatherSubField('headSwing', 'tailDir');
ota = odordilute.gatherSubField('headSwing', 'accepted');
ott = odordilute.gatherSubField('headSwing', 'maxTheta');
otv = odordilute.gatherField('ihtValid', 'headSwing', 'mean') == 1;
oteti = odordilute.gatherField('eti', 'headSwing', 'mean');


fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
left = (abs(otd) > 3*pi/4 & otv & oteti > 60);
right = (abs(otd) < pi/4 & otv & oteti > 60);
up = (abs(otd - pi/2) < pi/4 & otv & oteti > 60);
down = abs(otd + pi/2) < pi/4 & otv & oteti > 60;
bendpos = ott > pi/4;
bendneg = ott < -pi/4;
lp = mean(ota(left & bendpos));
ln = mean(ota(left & bendneg));
upos = mean(ota(up & bendpos));
un = mean(ota(up & bendneg));
rp = mean(ota(right & bendpos));
rn = mean(ota(right & bendneg));
dp = mean(ota(down & bendpos));
dn = mean(ota(down & bendneg));

bar([ln lp; un upos; rn rp; dn dp]);
set(gca, 'XTickLabel', {'headed to 0 degrees', 'headed to 90 degrees', 'headed to -180 degrees','headed to - 90 degrees'});
legend ('body angle < 0', 'body angle > 0'); 
title ('head swing accpetance odor');

fignum = fignum + 1; figureForPrinting(fignum); clf(fignum)
t = [odordilute.expt.track];
hs = [t.headSwing];
%size(hs)
rnhs = hs(right & bendneg);
hsnum = 0;
for j = 1:10
    for k = 1:10
        hsnum = mod(hsnum + 1,length(rnhs));
        illustrateHeadSwing(rnhs(hsnum), 'position', [40*j; 40*k], 'LineWidth', 2); hold on
    end
end
axis ([0 500 0 500]); axis equal