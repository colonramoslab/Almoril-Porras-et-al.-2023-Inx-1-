if (~exist('n2_cryo', 'var'))
    n2_cryo = ExperimentSet.fromMatFiles('E:\worm thermotaxis bin files\15C\N2\matfiles\n2_15c');
    alltracks = [n2_cryo.expt.track];
    so = alltracks(1).so;
    so.speedEndSTThresh = 0.0050;
    [alltracks.so] = deal(so);
    n2_cryo.executeTrackFunction('segmentTrack');
    allreo = [alltracks.reorientation];
    allst = [alltracks.sharpTurn];
end

   
if (~exist('tt', 'var'))
    tt = n2_cryo.gatherField('totalTime');
end

figure(1); clf();
alltracks(tt > 1200).plotPath('displacement', 'k-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5]);
hold on
rp = 8.9/2;
th = linspace(0,2*pi, 360);
plot (rp*cos(th), rp*sin(th), 'm--', 'LineWidth', 4);
trackinds = [118 158 239 248];
c = 'rgbk';
for j = 1:4
    alltracks(trackinds(j)).plotPath('displacement', c(j), 'LineWidth', 2);
end
axis equal

