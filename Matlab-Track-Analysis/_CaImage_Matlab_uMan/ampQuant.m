function [ampMat, fh] = ampQuant(dirs, incFrames, saveName, varargin)
%ampQuant Quantify response frequency and amplitude for files listed in
%dirs
%   Inputs:
%       dir, cell of strings pointing to file location
%       incFrames, array with values included in analysis, 
%             e.g. 1:10 for first 10 frames
%       saveName, output base for saving data & figure outputs

%% outputs
% ampMat, columns are samples % rows are responses
%  - Number of animals sampled, n, is column size
%  - Number of responses per animal are rows with values
%  - Values in psoitions are amplitude of response with respect to initial
%  signal

% saved outputs: 
% 1. savename.fig, figure illustrating response calls
% 2. savename-Var, variables including:
%     a. fluorescence signal, fluorMat
%     b. temperature stimulus, tempMat
%     c. response start times, respMat
%     d. response peak times, peakMat
%     e. response amplitudes, ampMat
%         - ampMat columns are samples % rows are responses

derThresh=0.01

varargin=assignApplicable(varargin);


%% First Extract response starts, respMat, and peaks, peakMat.
[tsT, tsF] = loadTS(dirs);
fluorMat=tsF.data(incFrames,:);
tempMat=tsT.data(incFrames,:);
tsF2=timeseries(tsF.data(incFrames,:),tsF.time(incFrames));
tsT2=timeseries(tsT.data(incFrames,:),tsT.time(incFrames));
[ ~, ~, respMat, peakMat]  = respFilter(tsF2, tsT2,1:length(incFrames),'derThresh',derThresh);


%% make a figure showing responses on heatmap, mainly QC
fh=figure(); imagesc(fluorMat')
hold on
% plot response starts
pMat=respMat;
for ii=1:size(pMat,2)
    tt=pMat(:,ii);
    tp=find(tt);
    for jj=1:length(tp)
        plot(tp(jj),ii,'linestyle','none','marker','o','markersize',5,'markerfacecolor','w');
%         for jj=1:length(stRep)
%             [~] = line('XData',stRep(jj),'YData',ii,...
%                 'Color','w','Marker','o', 'MarkerSize',10,'lineStyle','none')
    end
end
% plot peaks
peakPoints=peakMat>0;
pMat=peakPoints;
for ii=1:size(pMat,2)
    tt=pMat(:,ii);
    tp=find(tt);
    for jj=1:length(tp)
        plot(tp(jj),ii,'linestyle','none','marker','o','markersize',5,'markerfacecolor','g');
%         for jj=1:length(stRep)
%             [~] = line('XData',stRep(jj),'YData',ii,...
%                 'Color','w','Marker','o', 'MarkerSize',10,'lineStyle','none')
    end
end


% save plot
pName=strcat(saveName,'-ampPlot');
[ ~, ~ ] = vectorSave( fh, pName);
saveas( fh, pName,'fig');
close(fh);


%% Generate matrix of amplitudes, ampMat
% columns are samples
% rows are responses
% NaN is placefiller
% values indicate response amplitude relative to start of signal

% Define dimensions of ampMat
w=size(respMat,2); % Number of samples
h=max(sum(respMat)); % maximum number of responses for any sample
ampMat=NaN([w,h]);
try
for ii=1:w
    ff=fluorMat(:,ii);
    rr=respMat(:,ii);
    pp=peakMat(:,ii)>0;
    starts=find(rr);
    peaks=find(pp);
    for jj=1:length(starts)
        closePeaks=peaks-starts(jj);
        closestDiff=min(closePeaks(closePeaks>0));
        nextPeak=peaks(closePeaks==closestDiff);
        ampMat(ii,jj)=ff(nextPeak)-ff(starts(jj));
    end
end
catch ME
    ME
end
    
% save variables
vName=strcat(saveName,'-ampVars');
save(vName)

end

