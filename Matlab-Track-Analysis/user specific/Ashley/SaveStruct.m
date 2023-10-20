function [] = SaveStruct (structvar)
%function [] = SaveStruct (structvar)
%
%Saves the struct variable "structvar" to a .mat file that is located in
%the directory you choose using a dialog. Use LoadStruct to then return
%the variable to Matlab.

	
    %first go to the place where the file should be saved then save file
    %and then cd back
    currdir=cd;
    [FileName,PathName]=uiputfile('.mat', 'Create file to write');
    cd(PathName);
    save(FileName, 'structvar');
    cd(currdir);
    
    
end