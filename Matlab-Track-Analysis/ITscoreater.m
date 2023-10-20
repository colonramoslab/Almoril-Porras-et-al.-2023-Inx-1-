function [fPs,hRat,ITscore] = ITscoreater(dDir)
%ITater Call analyzeIT to examine isothermal tracking for contents of
%folder and subfolders
%   Detailed explanation goes here

% if inappropriate input, ask for file.
if ~exist('dDir','var')||isempty(dDir)
    [dDir]=uigetdir('C:\Users\jshha\Dropbox\ColonRamosLab\ITstory\ThermotaxisData');
end


% make file list from fPath
fL=dir([dDir, '/**/*.mat']);

% save to common 'parent' directory folder
saveDir=fullfile(dDir,'plots');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end
fPs={}; % store x-wise ratio of IT.
hRat={}; % store x-wise ratio of IT.
ITscore={}; % total fraction of time IT.
meanIT=[];

for ii=1:length(fL)
    fPath=fullfile(fL(ii).folder,fL(ii).name);
    [fPs{ii}, hRat{ii}, ITscore{ii}] = analyzeITscore(fPath,...
        'prefix',fL(ii).name(1:19),'saveDir',saveDir);
    meanIT(ii)=nanmean(ITscore{ii});
end




end