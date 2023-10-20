function [tsT,tsF, outDir] = retroPluckdStruct(fName, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Definitions
% Input
%   fName; % original data structure source location, with file name.

% optional inputs
% outName=''; % where to place the data, by default new folder with dStruct name in original folder
% ZOI=1; % column extracted from data structure, usually distinct zones of neurons.

% outputs
% tsF, double timeseries with fluorescence info, column is sample.
% tsT, double timeseries with temperature information
% outDir, if not defined by user, created to store variables.

%% initialize
% optionals inputs
outDir=''; % where to place the data, by default new folder with dStruct name in original folder
ZOI=1; % column extracted from data structure, usually distinct zones of neurons.
% outputs
tsF=timeseries;
tsT=timeseries;


% varargin=assignApplicable(varargin)

%% QC & define empty optionals
% check to see if fName exists as matlab variable
if ~exist(fName)
    error('retroPluckdStruct: file, fName, does not exist')
elseif ~(exist(fName)==2)
    error('retroPluckdStruct: file location, fName, is not formatted appropriately')
end

% define output dir if not specified
if isempty(outDir)
    [bD, bF, ~]=fileparts(fName);
    outDir=fullfile(bD,bF);
end

% make output dir if needed
if ~exist(outDir)
    mkdir(outDir)
end


%% Load datastructure
aStruct=load(fName);
aField=fieldnames(aStruct);
vName=aField{1};
aStruct=aStruct.(vName);

%% Extract
tsF=timeseries;
tsT=timeseries;
tsF.time=aStruct(1).green.time./1000;
tsT=aStruct(1).temperature;
tsT.time=tsT.time./1000;


for ii=1:length(aStruct)
    tD=aStruct(ii).green.data(:,ZOI);
    tD=(tD-min(tD))./min(tD); % store dF/F
    try
    if length(tsF.data)==0
        tsF.data(1:length(tD),ii)=tD;
    elseif length(tD)>length(tsF.data)
        tsF.data(:,ii)=tD(1:length(tsF.data));
    else
        tsF.data(1:length(tD),ii)=tD;
    end
    catch ME
        fy=1;
    end
end

% figure(); plot(tsF)
save(fullfile(outDir,'fluoStim.mat'),'tsF','tsT');
tsFout=[tsF.time, tsF.Data];
save(fullfile(outDir,'time_fluo.txt'),'tsFout','-ascii');
tsTout=[tsT.time, tsT.Data];
save(fullfile(outDir,'time_Temp.txt'),'tsTout','-ascii');



end

