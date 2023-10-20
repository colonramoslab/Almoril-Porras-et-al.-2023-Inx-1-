function [ sTimeSec, sTemp ] = getStimInfo(stimFile,stimDir)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Read timing information from stimulus file
% get stimulus file.
% MAKE AUTO BASED ON MIN DIFF FROM START TIME OF VIDEO



if isempty(stimFile)|isempty(stimDir)
    [stimFile,stimDir]=uigetfile('*.txt');
end

sInfo=readtable(fullfile(stimDir,stimFile));
sTime=nan([size(sInfo,1),3]);

% Extract start time contained in filename for stim file
sTime(1,1)=str2double(stimFile(8:9));
sTime(1,2)=str2double(stimFile(10:11));
sTime(1,3)=str2double(stimFile(12:13))+ str2double(stimFile(15:17))/1000;

sTime(:,3)=sInfo{:,4};

% Increment minutes based on seconds in column 4 of sInfo
testVar=sInfo{:,4};
incP=find(diff(testVar)<0);
incMat=zeros([size(sInfo,1),1]);
for i=1:length(incP)
    incMat(incP(i)+1:end)=incMat(incP(i)+1:end)+1;
end

sTime(:,2)=sTime(1,2)+incMat;

% Increment hours based on minutes in column 4 of sInfo
testVar=sTime(:,2);
incP=find(diff(testVar)<0);
incMat=zeros([size(sInfo,1),1]);
if ~isempty(incP)
    for i=1:length(incP)
        incMat(incP(i)+1:end)=incMat(incP(i)+1:end)+1;
    end
end

sTime(:,1)=sTime(1,1)+incMat;

% Convert to common sec-based frame
sTime(:,1)=sTime(:,1)*60*60;
sTime(:,2)=sTime(:,2)*60;
sTimeSec=sum(sTime,2);

sTemp=sInfo{:,1};
end



