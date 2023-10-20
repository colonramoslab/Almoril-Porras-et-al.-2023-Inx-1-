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
saveDir=fullfile(dDir,'IT_tracks_scores');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end
fPs={}; % store x-wise ratio of IT.
hRat={}; % store x-wise ratio of IT.
ITscore={}; % total fraction of time IT.
meanIT=[];

for ii=1:length(fL)
    fPath=fullfile(fL(ii).folder,fL(ii).name);
    expname=split(fL(ii).name,'_w1a_tracks');
    [fPs{ii}, hRat{ii}, ITscore{ii}] = analyzeITscore(fPath,...
        'prefix',expname(1),'saveDir',saveDir);
    meanIT(ii)=nanmean(ITscore{ii});
end




end