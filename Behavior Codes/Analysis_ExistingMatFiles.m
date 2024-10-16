dDir='U:\ThermotaxisGroup\BasicTTXanalysis\Tc15';
fileL=dir([dDir, '/*.mat']);


% Looking for 'preferred angle' that would have maximal travel before
% reorientation.

exptDir=cell(size(fileL)); % Store all direction of each run for expt as cell
exptDur=cell(size(fileL)); % Store all Duration  of each run for expt as cell
for ii=1:length(exptDur);
% Load & segment file
load(fullfile(fileL(ii).folder,fileL(ii).name));
segmentTracks(experiment_1)
% Confirm data. Skip 
[fignum]=ExpPlotTracks2(experiment_1);
% Extract direction & length for all *runs* of each *track*
trackDir=cell(size(experiment_1.track));
trackDur=cell(size(experiment_1.track));
for jj=1:length(trackDir)
    trackDir{jj}=zeros(size(experiment_1.track(jj).run));
    for kk=1:length(trackDir{jj})
    trackDir{jj}(kk)=experiment_1.track(jj).run(kk).meanTheta;
    trackDur{jj}(kk)=experiment_1.track(jj).run(kk).runTime;
    end
end
exptDir{ii}=combineGroups(trackDir); % each cell contains matrix run# x track#
exptDur{ii}=combineGroups(trackDur); % each cell contains matrix run# x track#

% Plot & save each.
figure(); polaraxes; polarplot(exptDir{ii},exptDur{ii},'.k')
title(fileL(ii).name);
saveas(gcf, fullfile(fileL(ii).folder,[fileL(ii).name,'polPlot.fig']))

figure(); 
plot(exptDir{ii},exptDur{ii},'.k')
title(fileL(ii).name);
saveas(gcf, fullfile(fileL(ii).folder,[fileL(ii).name,'linPlot.fig']))

% Discretized plot, dir by duration
N=10;
[Y,E] = discretize(exptDir{ii},0:pi/10:pi);
m=nan([N,1]);
mSTD=nan([N,1]);
mCnt=nan([N,1]);

for kk=1:length(E)
    m(kk)=nanmedian(exptDur{ii}(Y==kk));
    mSTD(kk)=nanstd(exptDur{ii}(Y==kk));
    mCnt(kk)=sum(sum(Y==kk));
end
xPt=cumsum(diff(E))-mean(E(1:2));
mE=1.96*mSTD./sqrt(mCnt);
figure(); hold on;
plot(xPt,m(1:end-1), '-b','LineWidth',5); 
plot(xPt,m(1:end-1)+mE(1:end-1), '-b','LineWidth',1); 
plot(xPt,m(1:end-1)-mE(1:end-1), '-b','LineWidth',1);

title(fileL(ii).name);
saveas(gcf, fullfile(fileL(ii).folder,[fileL(ii).name,'discPlot.fig']))

end

% extract from NaN-buffered arrays in cells... How many total runs?
runCnt=0;
for ii=1:length(exptDur)
    runCnt=runCnt+sum(sum(~isnan(exptDur{ii})));
end

runDir=nan([runCnt,1]);
runDur=nan([runCnt,1]);
fRun=0;
for ii=1:length(exptDur)

    runs=~isnan(exptDur{ii});
    runNum=sum(sum(runs));
    
    fRun=fRun+1;
    lRun=fRun+runNum-1;
    
    runDir(fRun:lRun)=exptDir{ii}(runs);
    runDur(fRun:lRun)=exptDur{ii}(runs);
    fRun=lRun;
end

figure(); histogram(runDir); ylabel('Direction of travel (rad)');
saveas(gcf,'AllRuns_dirHist.fig');
figure(); histogram(runDur); ylabel('Duration of travel (sec)');
saveas(gcf,'AllRuns_durHist.fig');
runDirAb=abs(runDir);
figure(); plot(runDirAb,runDur,'.k')
yLs=get(gca,'ylim'); line([pi/2,pi/2], yLs)
set(gca,'xlim',[0,pi],'xdir','reverse','xtick',xPt)
xlabel('Direction (rad)');
saveas(gcf,'AllRuns_dirByDur.fig');

figure(); 
polaraxes; hold on
polarplot(runDirAb,runDur,'.k');
polarplot([0,nanmean(runDirAb)],[0,nanmean(runDur)],'-r')
polarplot([0,mode(runDirAb)],[0,mode(runDur)],'-b')
polarplot([0,nanmedian(runDirAb)],[0,nanmedian(runDur)],'-g')
saveas(gcf,'AllRuns_abPolarPlot.fig');

figure(); 
polaraxes; hold on
% polarplot(runDir,runDur,'.k');
polarplot([0,nanmean(runDir)],[0,nanmean(runDur)],'-r')
polarplot([0,mode(runDir)],[0,mode(runDur)],'-b')
polarplot([0,nanmedian(runDir)],[0,nanmedian(runDur)],'-g')
saveas(gcf,'AllRuns_PolarPlot.fig');

N=10;
[Y,E] = discretize(runDirAb,0:pi/N:pi);
m=nan([N,1]);
mSTD=nan([N,1]);
mCnt=nan([N,1]);
for ii=1:length(E)
    m(ii)=nanmean(runDur(Y==ii));
    mSTD(ii)=nanstd(runDur(Y==ii));
    mCnt(ii)=sum(sum(Y==ii));
end
xPt=cumsum(diff(E))-mean(E(1:2));
mE=1.96*mSTD./sqrt(mCnt);
figure(); hold on;
plot(xPt,m(1:end-1), '-b','LineWidth',5); 
plot(xPt,m(1:end-1)+mE(1:end-1), '-b','LineWidth',1); 
plot(xPt,m(1:end-1)-mE(1:end-1), '-b','LineWidth',1);
yLs=get(gca,'ylim'); line([pi/2,pi/2], yLs)
set(gca,'xlim',[0,pi],'xdir','reverse','xtick',xPt)
xlabel('Direction (rad)');
saveas(gcf,'AllRuns_dirByDur_Binned.fig');

