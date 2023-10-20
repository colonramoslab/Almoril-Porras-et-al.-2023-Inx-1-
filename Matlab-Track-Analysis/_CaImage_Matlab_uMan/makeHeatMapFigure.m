function [ fM,clims,ah ] = makeHeatMapFigure( tsT, tsF, varargin)
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

fWidth=300; fHeight=900;

tempRange=[17.5, 22.5];
interval=60; % marker position spacing
orderOpt=1;
quantOpt=0;
saveOpt=1;
spMark=1; % 1 for marking responses, 0 for not.
derThresh=0.01; % Threshold derivative for response detection: 0.015 for AIY (190604), 0.002 for AFD (190604)

meanPanels=1;
tempPanels=1;
hmPanels=2;

semOpt=1; % use SEM as error in mean trace
tempColor=[0, 0, 0];
clims=[];
maxResp=1.4;

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

set(fM, 'position',[560,50,fWidth,fHeight]);

%% Panel A: Average trace
meanPlace=1:meanPanels;
ah{1}=subplot(vertPanels,1,meanPlace);
hold on;
incResp=nan([length(incFrames),size(tsF.data,2)]);
for ii=1:size(tsF.data,2)
    incResp(:,ii)=tsF.data(incFrames,ii);
end
meanResp=nanmean(incResp,2);
cntResp=size(incResp,2);
respError=1.96.*nanstd(incResp,1,2)./sqrt(cntResp);
ciResp=[meanResp+respError, meanResp-respError];
pt=plot(tsT.time(incFrames),meanResp,'color',tempColor,'linewidth',3);
plot(tsT.time(incFrames),ciResp,'color',tempColor,'linewidth',1);

axis('tight')
set(gca,'ylim',ylims,'xlim',[tsT.time(incFrames(1)),tsT.time(incFrames(end))]);
ylabel('mean \DeltaF/F');
xlabel('time (s)');


%% Panel B: temperature protocol
tempPlace=(meanPanels+1):(meanPanels+tempPanels);
ah{2}=subplot(vertPanels,1,tempPlace);
%plot(tsT.time(incFrames)-tsT.time(1),tsT.data(incFrames),'color',tempColor);
plot(tsT.time(incFrames),tsT.data(incFrames),'color',tempColor); % Changed to allow discontinuous displays of interval data
axis('tight')
ylabel('temperature (C)');
xlabel('time (s)');
t1=ah{1}.Position;
t2=ah{2}.Position;
set(ah{2},'Position',[t1(1),t2(2),t1(3),t1(4)]);


%% Panel C: heatmap
hmPlace = (tempPlace(end)+1):(tempPlace(end)+hmPanels);
ah{3}=subplot(vertPanels,1,hmPlace);
dataMat=tsF.data(incFrames,:)';
if orderOpt==1;
    [ dataMat ] = orderDataMat( dataMat ); % order responses based on MaxInt
elseif orderOpt==2; % max of middle two quartiles
    mat2=size(dataMat,2);
    rankMat=dataMat(:,mat2/4:3*mat2/4);
    matRank=transpose(max(transpose(rankMat)));
    [Y,I]=sort(matRank);
    dataMat=(dataMat(I,:));
end

dataMat=flipud(dataMat);
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
    [ ~, ~, respCnt]  = respFilter(tsF2, tsT2,1:length(incFrames),'derThresh',derThresh);
%     respMat=respCnt(incFrames,:);
    figure(fM); % return attention to heatmap
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

ylabel('animals', 'fontsize', 18);
set(gca,'FontSize', 16)
set(gca,'XTickLabel','');


ah2=gca;ah1=ah{1};
pI1=ah1.Position;
pI2=ah2.Position;
pI2(1)=pI1(1); pI2(3)=pI1(3);
set(ah2,'Position',pI2);

colormap jet;
set(gca,'FontName','arial')


end

