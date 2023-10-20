function [ fluorMat, tempMat, sPeaks, peakMat] = respFind( fluor, temp, varargin )
%respFind Create logical index into fluorescence values, fluor, where
%responses occur (single point or range?)
%   Detailed explanation goes here

% pass a screen to only select samples from a particular window
% see makeScreen.m for help.
screen=ones([size(fluor)]);
sampWin=50; % @10hz=5sec.
derThresh=0.01; % 0.015 for AIY (190604), 0.002 for AFD (190604
normOpt=0; % normalize fluorescence values prior to analysis.

varargin=assignApplicable(varargin);

%% Normalize fluorescence to reduce dependence on signal magnitude
if normOpt
    fluor=fluor./max(fluor);
end

%% using derivative
derF=diff(fluor);
% set threshold @ mean + 3*SD
if isempty(derThresh)
    derThresh=mean(nanmean(derF))+3*mean(nanstd(derF));
end
% find peaks based on derivative
passThresh=and(derF>derThresh,screen(2:end,:));

% remove any singleton points, where
% instantaneous passThresh values... d''=2
d2Mat=diff(diff(passThresh));
testVal=2;
d2Fail=abs(d2Mat)==testVal;
pad=zeros([1,size(d2Fail,2)]);
d2Fail=[pad;d2Fail;pad];
d2Fail=logical(d2Fail);
passThresh(d2Fail)=0;

% find start
sPeaks=diff(passThresh);
b=zeros([1,size(sPeaks,2)]);
sPeaks=[b;sPeaks;b]; % corrects for absence of terminal Diff
sPeaks=sPeaks>0; % converts to logical for start

% Add selection based on dF/F???
% sPeaks=findMax(fluor, sPeaks, derThresh);

% % QC1: plot with threshold
%  for i=1:size(passThresh,2)
%      figure(); plot(fluor(:,i));
%      yyaxis right; hold on; plot(passThresh(:,i)); plot(screen(:,i));
%  end
% No false positives in test wt data.
% A few false negatives (~5/40) with 4/5 in lower quality sample
% May revisit, but not important ATM.
% COULD lower threshold to 2*STD, then use peak width to improve quality

% % QC2: plot with threshold
% for i=1:size(passThresh,2)
%     figure(); plot(fluor(:,i)); yyaxis right; plot(sPeaks(:,i));
% end

% extract +/- 10sec from each peak
[fluorMat, tempMat, peakMat ] = pullPoints2( fluor, temp, sPeaks, sampWin );

% testPlots_PeaksStarts.m





% % QC4:
% mP=mean(respMat,2);
% figure(); plot(mP)



end

