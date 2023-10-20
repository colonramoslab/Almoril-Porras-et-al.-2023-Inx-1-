function [folderList,outFail, ME] = caImageSubfolders(parentFolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ME=struct;
% GETTING A DIRECTORY FOR LATER USE.
if isempty(parentFolder)
    parentFolder=uigetdir('Select folder to create list of experiments:');
end

cd(parentFolder);

% get all subdirectories.
pFcontents=dir(parentFolder); % pull all folders (&files) from parent
folderList=cell([length(pFcontents),1]); % initialize to max length possible
pFcnt=0; % counter of valid files
outFail=0; % how many folders failed to have a valid video, if fail, look inside later

% select only folders, fill serially based on counter pFcnt
for i=1:length(pFcontents)
    if exist(pFcontents(i).name,'dir')==7 % could use isfolder in later versions.
        if ~or(strcmp(pFcontents(i).name,'.'),strcmp(pFcontents(i).name,'..'))
            pFcnt=pFcnt+1;
            folderList{pFcnt}=fullfile(parentFolder,pFcontents(i).name);
        end
    end
end

folderList=folderList(1:pFcnt); % get rid of empty positions

% initialize info for comparisons
sampNum=length(folderList);
ListOLists=cell([sampNum,1]);

% Process videos within each folder
for i=1:sampNum
    try
        posNum=i;
        folderName=folderList{i};
        fileName=getImageFile(folderName);
        if ~isempty(fileName)
            [~, ~] = CaImageAnalysis(fileName);
        else
            outFail=outFail+1;
        end
    catch ME
        ME
    end
end



end

