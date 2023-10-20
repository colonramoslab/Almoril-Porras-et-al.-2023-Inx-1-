function [comboMat1,comboTime1,comboMat2,comboTime2] = responseCompare(dirs1,dirs2,legValues, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Necessary preceding analyses:

% dirs1: for group 1 list of data directories as cell array
% dirs2: for group 2 list of data directories as cell array
% legValues: cell array of group names, {'group1','group2'}

%% Initialize optional inputs.
% add path to functions to allow this to work: CaImaging_uMan folder
dirPathAdd='C:\Users\jshha\Dropbox\ColonRamosLab\Matlab\CaImaging_uMan\_CaImage_Matlab_uMan';

saveDir=pwd;

% Global Analysis parameters
%  Precursor to individua response analysis
%  Includes heatmap of responses with calls indicated
windO='Temp';   % Save name component for figures of data
frameRate=0.1;  % Can calculate as mean(diff(fluor.time)) below. 10Hz here.

incFrames={}; % overall window of the protocol (in frames), defaults below if not input
incFrames1=[1000:1900];
incFrames2=[1000:1900];  % allows specifying two windows, e.g. of same dirs.
% This means that global figures and all responses will be analyzed
% from 1000th frame to 1900th frame.

maxResp=.8; % Maximum on y-axis of average response trace
clims=[.05,1.3];    % Sets limits for colormap of heatmap
spM=1;  % Indicates that responses will be shown on heatmap
% posInfo=[5.5,1,3,5];

% Response parameters
sampWin=100; % sample window (frames or seconds?) before/after calcium response to include
cLrs={'k','r'};





varargin=assignApplicable(varargin)



if exist(dirPathAdd,'dir')
    addpath(genpath(dirPathAdd));   % Probably unnecessary, but add path to functions
end

if isempty(legValues)
    legValues={'group1','group2'};
end

if isempty(incFrames)
    incFrames{1}=incFrames1;
    incFrames{2}=incFrames2;
end

% cell array of matrices with response information
respMats={};
tempMats={};
respCnts={};
RRs={};
gCnt=0;

%
gCnt=1;
[temp, fluor] = loadTS(dirs1);
[ respMats{gCnt}, tempMats{gCnt}, respCnts{gCnt}, RRs{gCnt}] = analyze_FreqAmp(dirs1, incFrames{1},'sampWin',sampWin);

%
gCnt=2;
[temp, fluor] = loadTS(dirs2);
[ respMats{gCnt}, tempMats{gCnt}, respCnts{gCnt}, RRs{gCnt}] = analyze_FreqAmp(dirs2, incFrames{2},'sampWin',sampWin);

% Extract rise & fall of responses to eliminate noise from neighboring responses

% Sample 1.
% Prepare response & normalize by initial?
rM=respMats{1}; % For now, pick one response
rM=rM(size(rM,1)/2:end,:); % Cut out pre-response phase
% extract & separate clean rise & decay phases to allow fitting
[riseMat1,decayMat1,riseTime1,decayTime1] = riseDecaySplitter(rM);
% Put back together again for display purposes
[comboMat1,comboTime1] = phaseCombiner(riseMat1,riseTime1,decayMat1,decayTime1);

% Sample 2.
% Prepare response & normalize by initial?
rM=respMats{2}; % For now, pick one response
rM=rM(size(rM,1)/2:end,:); % Cut out pre-response phase
% extract & separate clean rise & decay phases to allow fitting
[riseMat2,decayMat2,riseTime2,decayTime2] = riseDecaySplitter(rM);
% Put back together again for display purposes
[comboMat2,comboTime2] = phaseCombiner(riseMat2,riseTime2,decayMat2,decayTime2);


%% PLOTS

% Plot raw data
figure(); hold on;
plot(respMats{1},['-',cLrs{1}])
plot(mean(respMats{1},2),['-',cLrs{1}],'linewidth',10)
plot(respMats{2},['-',cLrs{2}])
plot(mean(respMats{2},2),['-',cLrs{2}],'linewidth',10)
set(gca,'xlim', [0,2*sampWin])
legend([cLrs{1},'=',legValues{1}], [cLrs{2},'=',legValues{2}])
title('Calcium response profiles')
ylabel('Fluorescence \DeltaF/F');
xlabel('time (frames)');
saveas(gcf,fullfile(saveDir,['mixedResponseProfiles_',legValues{1},'-',legValues{2},'.fig']));
vectorSave(gcf,fullfile(saveDir,['mixedResponseProfiles_',legValues{1},'-',legValues{2},'.pdf']));

% plot outcome with comboMat
figure(); hold on;
% plot(riseTime15,riseMat15,'.k');
% plot(decayTime15,decayMat15,'--k');
plot(comboTime1,comboMat1,['-',cLrs{1}]);
plot(nanmean(comboTime1,2),nanmean(comboMat1,2),['-',cLrs{1}],'linewidth',10);
% plot(riseTime25,riseMat25,'.r');
% plot(decayTime25,decayMat25,'--r');
plot(comboTime2,comboMat2,['-',cLrs{2}]);
plot(nanmean(comboTime2,2),nanmean(comboMat2,2),['-',cLrs{2}],'linewidth',10);
set(gca,'xlim', [0,sampWin*frameRate])
legend([cLrs{1},'=',legValues{1}], [cLrs{2},'=',legValues{2}])
title('Individual calcium response profiles')
ylabel('Fluorescence \DeltaF/F');
xlabel('time (sec)');
saveas(gcf,fullfile(saveDir,['ResponseProfiles_',legValues{1},'-',legValues{2},'.fig']));
vectorSave(gcf,fullfile(saveDir,['ResponseProfiles_',legValues{1},'-',legValues{2},'.pdf']));


% Plot normalized to maximal response
normMat15=comboMat1./max(comboMat1,[],1);
normMat25=comboMat2./max(comboMat2,[],1);
figure(); hold on;
plot(comboTime1,normMat15,['-',cLrs{1}]);
plot(nanmean(comboTime1,2),nanmean(normMat15,2),['-',cLrs{1}],'linewidth',10);
plot(comboTime2,normMat25,['-',cLrs{2}]);
plot(nanmean(comboTime2,2),nanmean(normMat25,2),['-',cLrs{2}],'linewidth',10);
set(gca,'xlim', [0,sampWin*frameRate], 'ylim',[0,1.2])
legend([cLrs{1},'=',legValues{1}], [cLrs{2},'=',legValues{2}])
title('Normalized calcium response profiles')
ylabel('Fluorescence \DeltaF/F');
xlabel('time (sec)');
saveas(gcf,fullfile(saveDir,['normResponseProfiles_',legValues{1},'-',legValues{2},'.fig']));
vectorSave(gcf,fullfile(saveDir,['normResponseProfiles_',legValues{1},'-',legValues{2},'.pdf']));


end

