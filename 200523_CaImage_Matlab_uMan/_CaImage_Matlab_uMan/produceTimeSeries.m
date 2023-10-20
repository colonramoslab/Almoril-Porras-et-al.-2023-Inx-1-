function [ tsF, tsT ] = produceTimeSeries(imFile,imDir)
%produceTimeSeries Summary of this function goes here
%   Outputs:
%       tsF, timeseries of fluorescence
%       tsT, timeseries of temperature

% Extend to multiple types of results.csv files (e.g. for different zones
% or neurons)? Maybe setup loop through dir(fullfile(imDir,'*.csv')), then
% create tsF & tsT z-stack. Would need to handle in downstream
% operations... Alternately, could create separate tsF & tsT informed by
% results name...

%% Prior this function:
% Quantify ROIs from video file in imageJ
% 1. With first ROI as background, Create ROIs in ROI manager for all neurons
%    - analyze-> tools-> ROI manager
%    - save ROIs within ROI manager, more-> save
% 2. Set measurement conditions
%   - analyze-> set measurements-> mean gray value
% 3. Meaure with ROI measure, - more-> multimeasure
%   - measure all slices & one row per slice should be checked
% 4. Set save conditions within measurement table
%   - Results-> Options, uncheck all Results Table Options values
% 5. Save as default name, 'Results.csv', within same directory as video.

% insert guiSelect vid file and strip rest of the locations...


%% Assumes directory structure.
baseDir=fileparts(imDir);
stimDir=fullfile(baseDir,'stimFiles');
if ~exist(stimDir,'dir')
    baseDir=fileparts(fileparts(imDir));
    stimDir=fullfile(baseDir,'stimFiles');
end
if ~exist(stimDir,'dir')
    warning('You are missing the stimDir');
end



[vTimeSec,elapsedT,iDay] = getuManTime_Yale('imFull',fullfile(imDir,imFile));

d=strcat(iDay(6:7),iDay(9:10),iDay(3:4));
hT=vTimeSec(1)/(60*60);
h=floor(hT);
mT=(hT-h)*60;
m=floor(mT);
sT=(mT-m)*60;
s=floor(sT);

dStr=strcat(d,num2str(h),num2str(m),num2str(s));
vStr=str2num(dStr);
% goal. Find stim file with closest timing.
list=dir(stimDir);

stimFile='';
minDelta=100000;

% find closest stimFile in time to uMan start time (vStr)
% SOMETIMES FAILS, WHY?
for i=1:size(list,1)
    testStr=list(i).name;
    testVar=strsplit(testStr,'_');
    if str2num(testVar{1});
        numVar=str2num([testVar{1},testVar{2}]);
        testDelta=vStr-numVar;
        if abs(testDelta)<abs(minDelta)
            minDelta=testDelta;
            stimFile=testStr;
        end
    end
end

% Extract time & temp from stimFile

if isempty(stimFile)|isempty(stimDir)
    pTit=strcat('Please select stimulus for: ',imFile);
    [stimFile,stimDir]=uigetfile('*.txt',pTit );
end
[ sTimeSec, sTemp ] = getStimInfo(stimFile,stimDir);

% get results from ROIs, using mean gray value ATM with an initial
% 1st measurement is background reading for subtraction...
% save without column headers
if exist(fullfile(imDir,'Results.csv'))
    % should read first line... If text delete it, if first row is diff==1
    % just need to cycle open rather than csvread...
    fVals= csvread(fullfile(imDir,'Results.csv'));
else
    error('You must create a file Results.csv with fluorescence measurements');
end
if mean(diff(fVals(:,1)))==1 % saved with row values in file
    fVals=fVals(:,2:end);
end

%BGnormVal=fVals-repmat(fVals(:,1),[1,size(fVals,2)]);
% Can incorporate later.
F0=repmat(min(fVals),[size(fVals,1),1]);
dFtoF=(fVals-F0)./F0;



% create time series sampled at 10hz to ease alignment
tsF=timeseries(dFtoF,vTimeSec);
tsT=timeseries(sTemp,sTimeSec);
[tsF,tsT] = synchronize(tsF,tsT,'Uniform','Interval',.1);
save(fullfile(imDir,'fluoStim.mat'),'tsF','tsT');
tsFout=[tsF.time, tsF.Data];
save(fullfile(imDir,'time_fluo.txt'),'tsFout','-ascii');
tsTout=[tsT.time, tsT.Data];
save(fullfile(imDir,'time_Temp.txt'),'tsTout','-ascii');

end

