function [outFail] = figureWrapper(parentFolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nameRef='fluoStim.mat';
sName='HMwrap';

if isempty(parentFolder)
    parentFolder=uigetdir('Select folder to create list of experiments:');
end

% get all subdirectories.
cd(parentFolder)
pFcontents=dir(parentFolder); % pull all folders (&files) from parent

if ~exist('outFail','var')
    outFail=0; % how many folders failed
end

% select only folders, fill serially based on counter pFcnt
try
    for ii=1:length(pFcontents)
        if strcmp(pFcontents(ii).name,nameRef)
            load(nameRef, 'tsT', 'tsF');
            [ fM,~ ] = makeHeatMapFigure( tsT, tsF);
            saveName=fullfile(parentFolder,sName);
            [ ~, ~ ] = vectorSave( fM, saveName);
            saveas( fM, saveName,'fig');
            close(fM);
        elseif exist(pFcontents(ii).name,'dir')==7 % could use isfolder in later versions.
            if ~or(strcmp(pFcontents(ii).name,'.'),strcmp(pFcontents(ii).name,'..'))
                [outFail] = figureWrapper(fullfile(parentFolder,pFcontents(ii).name));
                cd(parentFolder)
            end
        end
    end
catch ME
    outFail=outFail+1;
end
end

