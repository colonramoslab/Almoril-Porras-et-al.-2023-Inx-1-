function [ temp_out, fluor_out ] = loadTS( dirs, varargin )
%loadTS Summary of this function goes here
%   Cell array of strings pointing to directories containing data to be
%   compiled

winRange=[];
bg=0; % 1= first recording is background sample. Added correction on 190820

varargin=assignApplicable(varargin);

TsMany=timeseries(); FsMany=timeseries();
for ii=1:length(dirs)
    if isfile(dirs{ii})
        [ta,tb]=fileparts(dirs{ii});
        tDir=fullfile(ta,tb);
        load(fullfile(tDir,'fluoStim.mat'));
    else
        load(fullfile(dirs{ii},'fluoStim.mat'));
    end
    if bg % background subtraction
        tsF.Data=tsF.Data(:,2:end)-tsF.Data(:,1);
    end
    TsMany(ii)=tsT;
    FsMany(ii)=tsF;
end


[ temp_out, fluor_out ] = combineTs( TsMany, FsMany );

if ~isempty(winRange)
    temp_out=getdatasamples(temp_out,winRange);
end

end

