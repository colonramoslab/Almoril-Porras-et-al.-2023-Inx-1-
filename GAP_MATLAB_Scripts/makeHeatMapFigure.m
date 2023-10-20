function [ fM,clims ] = makeHeatMapFigure( tsT, tsF, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% make figure with timeSeries list
% requires

% Add average trace at top

% Align starts of data better.


%% default values, can change by pairing with varargin
incFrames=[]; % vector specifying frames, e.g. 1:10 for first 10.
fM=[]; % old figure handle

sDir=pwd; % default to save in current working directory
fileName='heatmap'; % default file save name.
saveName=''; % will use sDir & filename if still empty later

fWidth=1000; fHeight=1600; %GAP default: fWidth=1600;

tempRange=[14.5 25.5];
interval=60; % marker position spacing
orderOpt=2;
quantOpt=0;
saveOpt=1;
spMark=0; % 1 for marking responses, 0 for not.

meanPanels=1;
tempPanels=1;
hmPanels=2;

semOpt=1; % use SEM as error in mean trace
tempColor=[0, 0, 0];
clims=[];
maxResp=1.5;

%% get any input values
varargin = assignApplicable (varargin);

ylims=[0,maxResp];

% if using cell formatting for multiple incFrame windows, unwrap them into
% single vector
if iscell(incFrames)
    if length(incFrames{1})==2 % format start and end only. Don't use but someone might by accident.
        tF=[];
        for ii=1:length(incFrames)
            tF=[tF,incFrames{ii}(1):incFrames{ii}(2)];
        end
        incFrames=tF;
    else % vector for each cell of incFrames
        tF=[];
        for ii=1:length(incFrames)
            tF=[tF,incFrames{ii}];
        end
        incFrames=tF;
    end
end
%% setup parameters
% new or old figure handle.
if isempty(fM)
    fM=figure();
end

% use all if not specified.
if isempty(incFrames);
    incFrames =1:numel(tsT.time);
end

% generate name from base values.
if isempty(saveName)
    saveName=fullfile(sDir,fileName);
end

vertPanels=meanPanels+tempPanels+hmPanels;

set(fM, 'position',[0,0,fWidth,fHeight]);%560,50

%% Panel A: Average trace
meanPlace=1:meanPanels;
ah1=subplot(vertPanels,1,meanPlace);
hold on;
incResp=nan([length(incFrames),size(tsF.data,2)]);
for ii=1:size(tsF.data,2)
    incResp(:,ii)=tsF.data(incFrames,ii);
end
meanResp=nanmean(incResp,2);
cntResp=size(incResp,2);
respError=1.96.*nanstd(incResp,1,2)./sqrt(cntResp);
ciResp=[meanResp+respError, meanResp-respError];
%for ii=1:size(tsF.data,2)
%    plot(tsT.time(incFrames),incResp(:,ii),'color',[0.83 0.83 0.83],'linewidth',0.5);
%end
plot(tsT.time(incFrames),meanResp,'color',tempColor,'linewidth',3);
plot(tsT.time(incFrames),ciResp,'color',tempColor,'linewidth',1);

axis('tight')
set(gca,'ylim',ylims,'Xticklabel',[]);
ylabel('GCaMP \DeltaF/F','fontsize', 16);
%xlabel('Time (s)');
%xticks(160:10:200);
%xticklabels({'0','10','20','30','40'})




%% Panel C: heatmap
hmPlace = (meanPanels+1):(meanPanels+hmPanels); %(tempPlace(end)+1):(tempPlace(end)+hmPanels); swapped line 107
subplot(vertPanels,1,hmPlace);
dataMat=tsF.data(incFrames,:)';
if orderOpt==1
    [ dataMat ] = orderDataMat( dataMat ); % order responses based on MaxInt
elseif orderOpt==2 % max of middle two quartiles
    mat2=size(dataMat,2);
    rankMat=dataMat(:,mat2/4:3*mat2/4);
    matRank=transpose(max(transpose(rankMat)));
    [Y,I]=sort(matRank);
    dataMat=(dataMat(I,:));
end

dataMat=flipud(dataMat);
% to remove BG:
%dataMat=dataMat(1:end-3,:);
if ~isempty(clims)
    im=imagesc(dataMat,clims);
else
    im=imagesc(dataMat);
    clims=get(gca,'CLim');
end

% mark response calls if applicable
% plot spike calls for each worm on heatmap
if spMark==1
    % make response calls, with re-ordered data
    tsF2=timeseries(dataMat',tsF.time(incFrames));
    tsT2=timeseries(tsT.data(incFrames,:),tsT.time(incFrames));
    [ ~, ~, respCnt]  = respFilter(tsF2, tsT2,1:length(incFrames));
%     respMat=respCnt(incFrames,:);
    for ii=1:size(respCnt,2)
        respList=respCnt(:,ii);
        stRep=find(respList);
        for jj=1:length(stRep)
            [~] = line('XData',stRep(jj),'YData',ii,...
                'Color','w','Marker','o', 'MarkerSize',10,'lineStyle','none')
        end
    end
end
%

axis('tight');
%     title(tiInfo);

ylabel('animals (n= )', 'fontsize', 16);
set(gca,'FontSize', 16)
set(gca,'XTickLabel','');

colorbar('EastOutside','XTickLabel',{'','','100%','','200%','','300%'},'XTick', 0:0.5:3);

ah2=gca;
pI1=ah1.Position;
pI2=ah2.Position;
pI2(1)=pI1(1); pI2(3)=pI1(3);
set(ah2,'Position',pI2);

colormap jet;
set(gca,'FontName','arial')

%% Panel B: temperature protocol
tempPlace=(hmPlace(end)+1):(hmPlace(end)+tempPanels); % (meanPanels+1):(meanPanels+tempPanels); swapped line 119
subplot(vertPanels,1,tempPlace);
%plot(tsT.time(incFrames)-tsT.time(1),tsT.data(incFrames),'color',tempColor);
plot(tsT.time(incFrames),tsT.data(incFrames),'color',tempColor); % Changed to allow discontinuous displays of interval data
axis('tight')
set(gca,'ylim');%,tempRange
ylabel('(\circC)','fontsize', 16);
xlabel('time (s)');
xticks(160:10:200);
xticklabels({'10','20','30','40','50'})
yticks(18:1:22);
yticklabels({'18',' ',' ',' ','22'})

end

