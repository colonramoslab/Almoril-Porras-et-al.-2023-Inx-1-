function [ fluorMat, tempMat, respCnt, RR, respMat] = analyze_FreqAmp(inArg, incFrames, varargin)
%analyze_FreqAmp Summary of this function goes here
%   Detailed explanation goes here

% default inArg is a cell array of string pointers to data
%   Alternatively, [temp,fluor] timeseries object valid

%incFrames= {1560:2160}; % could write routine to select from full figure
% would set to empty, then used data to make 1:end...
% should use cell array if creating multiple

%% optional inputs
sampWin=50; % Should use seconds here
derThresh=0.015; % 0.015 for AIY (190604), 0.002 for AFD (190604)
respIll=1; % QC step, show where responses are on heatmap

varargin=assignApplicable(varargin);

plPath='C:\Users\jshha\OneDrive\Documents\ColonRamosLab\140514_CaImaging\HelperFunctions\PlottingFunctions';
if exist(plPath,'dir')
    addpath(plPath);
end

if strcmp(class(inArg),'timeseries')
    temp=inArg(1); fluor=inArg(2);
else
    [temp, fluor] = loadTS(inArg);
end


%% re-format incFrames to index samples
frameSynth=[];
if iscell(incFrames)
    for i=1:length(incFrames)
        frameSynth=[frameSynth,incFrames{i}];
    end
elseif isnumeric(incFrames)
    if length(incFrames)==2
        frameSynth=incFrames(1):incFrames(2);
    else
        frameSynth=incFrames;
        incFrames={};
        incFrames{1}=[frameSynth(1),frameSynth(end)];
    end
else
    warning('Frames defined in incFrames are in the wrong format');
end

%% extract fluorescence within specified window
if incFrames{1}(2)>size(fluor.Data,1)
    incFrames{1}(2)=size(fluor.Data,1);
end

respMat=fluor.Data(incFrames{1},:);

%% find peaks using current derThresh & show on heatmap
if respIll % QC step, show where responses are on heatmap
    [ ~ ] = makeHeatMapFigure( temp, fluor,'incFrames',frameSynth,'spMark',1,'derThresh',derThresh);
end

%%  Calculate response frequency during relevant epoch
[ fluorMat, tempMat, respCnt]  = respFilter(fluor, temp, frameSynth,'derThresh',derThresh,'sampWin',sampWin);

responseCount=sum(respCnt);
sampDuration=0;
for i=1:numel(incFrames)
    W1=incFrames{i};
    sampDuration=sampDuration+ fluor.Time(W1(end))- fluor.Time(W1(1));
end
RR=responseCount./sampDuration;

