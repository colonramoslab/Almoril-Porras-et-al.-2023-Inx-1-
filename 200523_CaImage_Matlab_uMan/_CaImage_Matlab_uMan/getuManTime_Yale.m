function [vTimeSec,elapsedT,iDay] = getuManTime_Yale(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Get timing information from video file

imFull='';
approach='';

varargin=assignApplicable(varargin);
% select/get directory
if isempty(imFull)
    [imFile,imDir] =uigetfile({'*.tif'},'Select tif file with associated metadata');
    imFull=fullfile(imDir,imFile);
else
    [imDir,imFile,~] = fileparts(imFull);
end
nBase=strsplit(imFile,'.');


%% SWITCH/CASE BASED ON separate or embedded imageJ metadata or metamorph
% get names of image files in dir
%fList=dir(fullfile(imDir,'*.tif'));

% read metadata & output frame time in columns: hr, min, sec
H=imfinfo(imFull);
mDname=fullfile(imDir,strcat(nBase{1},'_metadata.txt'));
% Define approach based on available information.
if exist(mDname)
    approach='exMeta';
elseif isfield(H,'Software') %testmetamorph maybe with imfinfo
    approach='metaMorph';
else
    approach='inMeta';
end


switch approach
    case 'exMeta'
        tR1=fileread(mDname);
        %"Time": "2018-10-09 16:23:53 -0400",
        tR2=strsplit(tR1,'{');
        tR2=tR2(4:end);
        try
            [vTimeSec,elapsedT,iDay] = unpackvTime(tR2);
        catch % seems some changes to metadata order, but imfExtract may still work
            % may be able to do away with exMeta fork
            [vTimeSec,elapsedT,iDay] = imfExtract(imFull);
        end
    case 'inMeta'
        %% try imfinfo
        [vTimeSec,elapsedT,iDay] = imfExtract(imFull);
    case 'metaMorph'
        [vTimeSec,elapsedT,iDay] = mmExtract(imFull);
end



end

