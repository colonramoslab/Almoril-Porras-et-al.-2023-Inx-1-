function [ TsSingle, FsSingle ] = combineTs( TsMany, FsMany )
%combineTs Takes multiple fluorescence and temperature time series and
%combines using timeseries interpolation and making temperature unique to
%each fluorescence trace.
%   Detailed explanation goes here


% option 1. pull out each, synchronize to existing body, incorporate into
% body?

% option 2. pull out all. define important parameters. incorporate in
% ordered, but still pair-wise, manner. Prefered if some parameters
% influence best order or if better to only synchronize against single
% common reference.

%% option 1. Combine into single growing structure by making time relative.
TsSingle=timeseries; TsNew=timeseries; 
FsSingle=timeseries; FsNew=timeseries;

% Insert values from first data set. Normalize time to start
FsNew=FsMany(1);
FsSingle.data=FsNew.data; 
FsSingle.time=FsNew.time-FsNew.time(1);

% Insert values from first data set. Normalize time to start
% and Replicate temperature values. one for each fluorescence data point
TsNew=TsMany(1);
TsNew.data=repmat(TsNew.data,[1,size(FsNew.data,2)]);
TsSingle.data=TsNew.data; 
TsSingle.time=TsNew.time-TsNew.time(1);

for i=2:length(FsMany)
    % pull out next value set.
    FsNew=FsMany(i);
    TsNew=TsMany(i);
    TsNew.data=repmat(TsNew.data,[1,size(FsNew.data,2)]);
    % make time relative to start.
    TsNew.time=TsNew.time-TsNew.time(1);
    FsNew.time=FsNew.time-FsNew.time(1);
    % synchronize new data with existing data
    [TsSingle,TsNew]= synchronize(TsSingle,TsNew,'Uniform','Interval',.1);
    [FsSingle,FsNew]= synchronize(FsSingle,FsNew,'Uniform','Interval',.1);
    % incorporate new data into growing structure
    startC=size(TsSingle.data,2)+1;
    endC=size(TsSingle.data,2)+size(TsNew.data,2);
    TsSingle.data(:,startC:endC)=TsNew.data;
    FsSingle.data(:,startC:endC)=FsNew.data;
end
    
    


end

