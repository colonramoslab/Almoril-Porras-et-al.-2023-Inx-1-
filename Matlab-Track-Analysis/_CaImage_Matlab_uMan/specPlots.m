function [fh] = specPlots(respCnts, incFrames, tempProf,specWin,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% % Require inputs, 
% respCnts
% incFrames
% tempProf
% [x1,x2], specificity window 

gName='';
windO='';
samps=1:length(respCnts);
lStyle=repmat({'-'},[1,length(samps)]);
cCh={'b','k','r','g'};

varargin=assignApplicable(varargin);

while length(cCh)<length(samps)
    cCh=[cCh,cCh];
end

x1=specWin(1); x2=specWin(2);


% KSdensity: plot window showing stimulus & shading +/- 2*STD for responses
[ dMat] = combineGroupsFlip(respCnts);
fh=figure(); plot(tempProf);
ylabel('Temperature')
set(gca,'ylim',[16,21],'xlim',[incFrames(1)/10, incFrames(end)/10]);
set(gca,'fontname','arial','fontsize',20)
[~] = xTzero(fh)
% x1=nanmean(sS)-2*nanstd(sS); % from AFD
% x2=nanmean(sS)+2*nanstd(sS); % from AFD
y2=get(gca,'ylim');
hold on;
h1 = fill([x1/10 x1/10 x2/10 x2/10], [y2 fliplr(y2)], 'b','EdgeColor','none');
set(h1,'facealpha',.1)
  
for ii=1:length(samps)
    [ dMat] = combineGroupsFlip(respCnts(samps(ii)));
    respTime=findResp(dMat);
    rTs{ii}=respTime;
    yyaxis right;
    [f,xi] = ksdensity(respTime,'bandwidth',50);
    plot(xi/10,f,'linestyle',lStyle{ii},'color',cCh{ii},'marker','none')
end
title('AIY response kernal density')
set(gca,'xlim',[incFrames(1)/10, incFrames(end)/10]);
yyaxis right;
ylimR=get(gca,'ylim');
saveName=strcat('kDens_AIY-',windO);saveName=strcat(gName,'kDens',windO);
[ ~, ~ ] = vectorSave( gcf, saveName);
saveas( gcf, saveName,'fig');

% boxplot & plotEach
[ gArray ] = combineGroupsFlip( rTs);
saveName=strcat(gName,'boxPlot',windO);
figure(); boxplot(gArray);
set(gca,'ylim',[incFrames(1),incFrames(end)]);
xL=get(gca,'xlim'); yL=get(gca,'ylim');
view(90,90)
[ ~, ~ ] = vectorSave( gcf, saveName);
saveas( gcf, saveName,'fig');

saveName=strcat(gName,'plotEach',windO);
plotEachPointMean(gArray,'pntColors',{'b','k','r'});
set(gca,'ylim',yL);
view(90,90)
[ ~, ~ ] = vectorSave( gcf, saveName);
saveas( gcf, saveName,'fig');

% Penetrance, AIY Wt, two windows
saveName=strcat(gName,windO,'Penetrance-AFDlocked');
% plot penetrance within specified window, AFD-locked
selFrames=[floor(x1):ceil(x2)]; % AFD-locked window
RRs=cell([length(samps),1]); rT=respCnts(samps);
for ii=1:length(rT)
    RRs{ii}=rrDeRespCnt(rT{ii}, 'windO', selFrames)
end
[ dMat] = combineGroupsFlip(RRs);
[ fh ] = plotPenetrance(dMat,'compSets',[ones([length(samps)-1,1]),samps(2:end)']);
set(gca,'ylim',[0,1.2]);
[ ~, ~ ] = vectorSave( fh, saveName);
saveas( fh, saveName,'fig');


% plot penetrance within specified window, non-AFD-locked
saveName=strcat(gName,windO,'Penetrance-AFDlocked');
selFrames=[1000:floor(x1),ceil(x2):1900]; % non-AFD-locked window
RRs=cell([length(samps),1]); rT=respCnts(samps);
for ii=1:length(rT)
    RRs{ii}=rrDeRespCnt(rT{ii}, 'windO', selFrames)
end
[ dMat] = combineGroupsFlip(RRs);
[ fh ] = plotPenetrance(dMat,'compSets',[ones([length(samps)-1,1]),samps(2:end)']);
set(gca,'ylim',[0,1])
set(gca,'ylim',[0,1.2]);
[ ~, ~ ] = vectorSave( fh, saveName);
saveas( fh, saveName,'fig');

% TUNING METRIC. Per worm
% Fraction of responses within AFD window...
% VERY simple. Sum within X1 to X2
for ii=1:length(respCnts)
    t1=respCnts{ii};
    allR=sum(t1(incFrames,:),1);
    Ar=sum(t1(floor(x1):ceil(x2),:),1);
    AI{ii}=Ar./allR;
end
[ gArray ] = combineGroups( AI);

saveName=strcat(gName,'AFDindexEach',windO);
plotEachPointBar(gArray(:,samps),'pntColors',{'b','k','r'});
set(gca,'ylim',[0,1.2],'fontsize',20);
[ ~, ~ ] = vectorSave( gcf, saveName);
saveas( gcf, saveName,'fig');
[p,tbl,stats] = kruskalwallis(gArray(:,samps))
c = multcompare(stats)



end

