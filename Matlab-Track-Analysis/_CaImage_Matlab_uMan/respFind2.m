function [ fluorMat, tempMat, sPeaks, peakMat]= respFind2(fluor, temp, varargin )
%UNTITLED Summary of this function goes here
%   Developed 200420 b/c of noisy AFD data...


% Optional inputs
stMulti=1; % # of std above mean smoothed derivative. Good results with only 1.
figOpt=0;

screen=ones([size(fluor)]); % pass a screen to only select samples from a particular window
% see makeScreen.m for help.
sampWin=50; % @10hz=5sec.

varargin=assignApplicable(varargin);

%% Considering spike detection
fluor=fluor./max(fluor);

% movMean smoothing
fdSm=smoothdata(fluor,'movMean',20);

% take derivative
dfdSm=diff(fdSm);

% smooth derivative, improves performance
dfdSmSm=smoothdata(dfdSm,'movMean',20);

% create a derivative threshold based on distribution
derThresh=mean(dfdSmSm(:))+stMulti*std(dfdSmSm(:));
% Based on histogram below.

% Find response windows using this threshold
passThresh=dfdSmSm>derThresh;

% find start of response windows
sPeaks=diff(passThresh);
b=zeros([1,size(sPeaks,2)]);
sPeaks=[b;sPeaks;b]; % corrects for absence of terminal Diff
sPeaks=sPeaks>0; % converts to logical for start



%% Figures for QC
if figOpt
    figure(); plot(fdSm)
    title('movMean')
    
    figure(); plot(dfdSm)
    title('movMean dF')
    
    
    histogram(dfdSmSm); hold on;
    yL=get(gca,'ylim');
    line([derThresh,derThresh], yL)
    
    figure(); plot(dfdSmSm)
    xL=get(gca,'xlim');
    line(xL,[derThresh,derThresh])
    title('movMean dF smooth')
    
    
    detectRate=sum(logical(sum(passThresh,1)))/length(logical(sum(passThresh,1)))
    
    figure(); imagesc(fluor')
    figure(); imagesc(passThresh')
    figure(); imagesc(sPeaks')
    
end


% extract +/- 10sec from each peak
[fluorMat, tempMat, peakMat ] = pullPoints2( fluor, temp, sPeaks, sampWin );

end

