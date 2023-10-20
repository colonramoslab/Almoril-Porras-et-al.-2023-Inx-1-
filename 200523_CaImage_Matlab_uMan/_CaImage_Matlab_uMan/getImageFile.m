function fileName=getImageFile(parentFolder);
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fN='';
fileName='';
fType='.tif';

if isempty(parentFolder)
    parentFolder=uigetdir('Select folder to create list of experiments:');
end

% get files/dir in dir.
pFcontents=ls(parentFolder); % pull all folders (&files) from parent

% select image file
for i=1:size(pFcontents,1)
    tst=strfind(pFcontents(i,:),fType);
    if tst
        fN=pFcontents(i,1:tst+length(fType)-1);
        fileName=fullfile(parentFolder,fN);
        break;
    end
end


end

