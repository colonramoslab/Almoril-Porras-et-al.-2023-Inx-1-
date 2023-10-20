function [riseMat,decayMat,riseTime,decayTime] = riseDecaySplitter(rM,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


dThresh=0.0005; % minimal derivative value to be consider a 'rise' in calcium
frameRate=0.1;  % Can calculate as mean(diff(fluor.time)) prior to generating rM
% see 'Analysis_ResponseKinetics.m'

varargin=assignApplicable(varargin);



[len,samps]=size(rM); % will use these values in loop
riseTime=nan([len,samps]);  % Convert to seconds for fitting
decayTime=nan([len,samps]);  % Convert to seconds for fitting
riseMat=nan([len,samps]);   % place calcium rise phases here
decayMat=nan([len,samps]);   % place calcium decay phases here
tims=frameRate:frameRate:len*frameRate;

% Extract rising & decaying phases
%  Rise: until first derivative below threshold, lTh
%  Decay: until first derivative above threshold, hTh
lTh=0.0005;
hTh=lTh; % single threshold probably fine. considering

for ii=1:size(rM,2)
    
    % extract data to examine
    inst=rM(:,ii); % pull out single sample
    dInst=diff(inst);   % find rate of change for identifying phases.
    % Note: smoothing may be necessary depending on data quality.
    dInst=[dInst(1);dInst]; % assume first step has same rate as second
    
    % delimit phases
    % First rise is between indices riseStart & riseStop
    % First decay is between indices decayStart & decayStop
    rThreshd=(dInst>dThresh); % pass derivative threshold to be rising.
    rises=find(rThreshd);
    riseStart=rises(1);
    fThreshd=~rThreshd;
    falls=find(fThreshd);
    riseStop= falls(falls>riseStart);
    if isempty(riseStop)
        riseStop=length(inst);
    else
        riseStop=riseStop(1);
    end
    
    decayStart=riseStop;
    decayStop=rises(rises>decayStart);
    if ~isempty(decayStop)
        decayStop=decayStop(1);
    else
        decayStop=length(inst);
    end
    
    % assign values
    riseInds=riseStart:riseStop;
    riseTime(1:length(riseInds),ii)=tims(riseInds);
    riseMat(1:length(riseInds),ii)=inst(riseInds);
    decayInds=decayStart:decayStop;
    decayTime(1:length(decayInds),ii)=tims(decayInds);
    decayMat(1:length(decayInds),ii)=inst(decayInds);
    
end

end

