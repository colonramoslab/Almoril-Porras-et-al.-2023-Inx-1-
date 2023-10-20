function [fh, ll, lr] = respDensity(respCnts, samps, temp, incFrames,xWin)

% respCnts, array of cells with binary response marking
% samps, selection of samples to be taken from respCnt
% temp, linear array of temperature stimulus values, 10hz? 
% incFrames, sampling window for data

windO='ramp';
saveName=strcat('kDens_AIY-',windO);
lStyle=repmat({'-'},[length(samps),1]);
cCh={};
for ii=1:length(samps)
    cCh{ii}=[rand,rand,rand];
end

x1=xWin(1);
x2=xWin(2);

% WT, AIY: plot window showing stimulus & shading +/- 2*STD for responses
% [ dMat] = combineGroupsFlip(respCnts);
fh=figure(); ll= plot(temp);
ylabel('Temperature')
set(gca,'ylim',[16,21],'xlim',[incFrames(1)/10, incFrames(end)/10]);
set(gca,'fontname','arial','fontsize',20)
[~] = xTzero(fh);
% x1=nanmean(sS)-2*nanstd(sS); % from AFD
% x2=nanmean(sS)+2*nanstd(sS); % from AFD
y2=get(gca,'ylim');
hold on;
h1 = fill([x1/10 x1/10 x2/10 x2/10], [y2 fliplr(y2)], 'b','EdgeColor','none');
set(h1,'facealpha',.1)
    
for ii=1:length(samps)
    [ dMat] = combineGroupsFlip(respCnts(samps(ii)));
    respTime=findResp(dMat);
    yyaxis right;
    [f,xi] = ksdensity(respTime,'bandwidth',50);
    lr= plot(xi/10,f,'linestyle',lStyle{ii},'color',cCh{ii},'marker','none');
end
title('AIY response kernal density')
set(gca,'xlim',[incFrames(1)/10, incFrames(end)/10]);
yyaxis right;
ylimR=get(gca,'ylim');
[ ~, ~ ] = vectorSave( gcf, saveName);
saveas( gcf, saveName,'fig');
end

