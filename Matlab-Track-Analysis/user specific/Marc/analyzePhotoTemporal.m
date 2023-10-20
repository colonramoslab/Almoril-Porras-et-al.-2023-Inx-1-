basedir = '\\labnas1\Share\Ashley Extracted\temporal\';
nums60 = [1750, 1830, 1915, 1950];%, 1850, 2015];

for j = 1:length(nums60)
    fn60{j} = [basedir '20_' num2str(nums60(j)) '_tracks.bin']; %#ok<SAGROW>
end
if (~exist('photo60', 'var'))
    minpts = 500;
    photo60 = ExperimentSet.fromFiles(fn60{:},'minpts',minpts);
end
existsAndDefault('discardLast10mins', false);
if (discardLast10mins)
    photo60.executeExperimentFunction('trimTracks', -1:600,[]);
    discardLast10mins = false;
end
photo60.defaultTitle = 'Temporal Phototaxis 60 sec cycle';

speedthresh = 2.5;
for j = 1:length(photo60.expt)
    ms = photo60.expt(j).gatherField('speed','mean');
    photo60.expt(j).track = photo60.expt(j).track(ms > speedthresh);
end

for j = 1:length(photo60.expt)
    period = 240;
    f = 1:length(photo60.expt(j).elapsedTime);
    lighton = mod(f, 2*period) < period;
    lightoff = ~lighton;
    interval = 0.25;
    %{
    timeon = interval * mod(f,2*period);
    timeon(lightoff) = 0;
    timeoff = interval * (mod(f,2*period) - period);
    timeoff(lighton) = 0;
    %break;
    %}
    lighttime = interval*(mod((f+period),2*period) - period);
    %lighttime is the time since the light turned on (in the light)
    %and (minus) the time until the light turns on (in the dark)
    photo60.expt(j).executeTrackFunction('addGlobalQuantity', 'lighton', photo60.expt(j).elapsedTime, lighton);
    photo60.expt(j).executeTrackFunction('addGlobalQuantity', 'lightoff', photo60.expt(j).elapsedTime, lightoff);
    photo60.expt(j).executeTrackFunction('addGlobalQuantity', 'lighttime', photo60.expt(j).elapsedTime, lighttime);
end

fignum = 0;
existsAndDefault('makeSpeedFigs', true);

if (makeSpeedFigs)
    fignum = fignum + 1; figure(fignum); clf(fignum);
    photo60.meanField2vsField1('lighttime', 'speed', -60:2:60); hold on;
    %yl = get(gca, 'Ylim'); yl(1) = 0;
    yl = [0 5];
    plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);

    fignum = fignum + 1; figure(fignum); clf(fignum);
    photo60.meanField2vsField1('eti', 'speed', 0:10:1200);

    earlyinds = find(photo60.gatherField('eti') < 600);
    lateinds = find(photo60.gatherField('eti') > 600);

    fignum = fignum + 1; figure(fignum); clf(fignum);
    photo60.meanField2vsField1('lighttime', 'speed', -60:2:60, 'inds', earlyinds); hold on;
    %yl = get(gca, 'Ylim'); %yl(1) = 0;
    plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);
    title ('first 10 minutes');

    if (~isempty(lateinds)) 
        fignum = fignum + 1; figure(fignum); clf(fignum);
        photo60.meanField2vsField1('lighttime', 'speed', -60:2:60, 'inds', lateinds); hold on;
        plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);
        title ('last 10 minutes');
    end
end

badratio = photo60.evaluateTrackExpression('mean(track.getDerivedQuantity(''ihtValid''))');
fignum = fignum + 1; figure(fignum); clf(fignum);
hist(badratio); title ('histogram of fraction of track points that have invalid head tail');

existsAndDefault('fixExpt', false);
if (fixExpt)
    badlimit = 0.975;
    disp(['ok, fixing ' num2str(sum(badratio < badlimit)) ' bad tracks']);
    mra = MaggotReAnalyzer;
    mra.contourScale = 3;
    mra.fixExperimentSet(photo60, 'badlimit', badlimit);
    fixExpt = false;
end
    
existsAndDefault('fixHT', true);
if (fixHT)
    disp('fixing HT orientation');
    photo60.executeTrackFunction('fixHTOrientation','mintime',5);
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
fignum = fignum + 1; figure(fignum); clf(fignum);
photo60.makeReorientationHistogram('lighttime', -60:2:60);

fignum = fignum + 1; figure(fignum); clf(fignum);
binsize = 10;
tx = (-60 + binsize/2):binsize:(60-binsize/2);
hst = photo60.makeHistogram('lighttime', tx, 'hsstart');
allt = photo60.makeHistogram('lighttime', tx);
plot (tx, hst./allt * 240); xlabel ('light time'); ylabel ('rate (min^{-1})'); title ('headsweep initiation rate vs. cycle');
embiggen();

fignum = fignum + 1; figure(fignum); clf(fignum);
r = photo60.gatherField('lighttime', 'reo', 'start');
r2 = photo60.gatherField('lighttime', 'reo', 'end');
r2(r2 < r) = r2(r2 < r) + 120;
valid = (r2 - r) < 30;
rc = (r+r2)/2;
rc = mod(rc + 60,120)-60;
n = photo60.gatherSubField('reorientation', 'numHS');
[x,meany] = meanyvsx(rc(valid), n(valid), -60:10:60); plot (x,meany,'LineWidth',2); title ('mean num headsweeps / reorientation');
xlabel ('reorientation start time'); ylabel ('<num hs>'); embiggen

fignum = fignum + 1; figure(fignum); clf(fignum);
sd = photo60.gatherSubField('reorientation', 'prevDir');
ed = photo60.gatherSubField('reorientation', 'nextDir');
dt = diff(unwrap([sd;ed])); 
[x,meany] = meanyvsx(rc(valid), rad2deg(abs(dt(valid))), -60:10:60);
[x2,meany2] = meanyvsx(rc(valid & n > 0), rad2deg(abs(dt(valid & n>0))), -60:10:60);
[x3,meany3] = meanyvsx(rc(valid & n == 1), rad2deg(abs(dt(valid & n == 1))), -60:10:60);
plot (x,meany,'b-',x2,meany2,'r--',x3,meany3,'g-.','LineWidth',2); title ('mean |angle change| / reorientation'); ylim([0 120])
legend('all reorientations', 'reorientations with at least one headswing', 'reorientations with exactly one headswing','Location','Best');
embiggen();

fignum = fignum+1; figure(fignum); clf(fignum);
tx = [-180:10:180];
dtd = rad2deg(dt);
plot (tx, hist(dtd(valid & n == 0),tx)/sum(valid & n==0), 'b-', tx, hist(dtd(valid & n == 1),tx)/sum(valid & n==1), 'r--', tx, hist(dtd(valid & n>1),tx)/sum(valid & n>1), 'g-.','LineWidth',2);
title ('Distribution of Angle Change during a reorientation'); ylabel ('fraction');
legend ('0 headswings', '1 headswing', '>1 headswing');
embiggen();

fignum = fignum+1; figure(fignum); clf(fignum);
photo60.meanField2vsField1('lighttime', 'speed', -60:2:60, 'runs'); hold on;
plot ([0 0], yl, 'r--','LineWidth',2); ylim(yl);
title ('speed in runs only'); hold off

fignum = fignum+1; figure(fignum); clf(fignum);
[x,meany] = meanyvsx(rc(valid), r2(valid)-r(valid), -60:10:60);plot (x,meany,'LineWidth',2); title ('mean reorientation duration');
xlabel ('reorientation start time'); ylabel ('<duration (s)>'); embiggen

fignum = fignum+1; figure(fignum); clf(fignum);
photo60.makeHeadSwingAcceptanceHistogram('lighttime', -55:10:55, 'atend', true, 'LineWidth',2); embiggen();