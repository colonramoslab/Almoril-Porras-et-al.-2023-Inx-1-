function [structvar] = LoadStruct ()
%function [structvar] = LoadStruct ()
%
%Loads the struct variable you specify in a dialog from a .mat file.

	

    
    
    %------load the struct
    currdir=cd;
    [FileName,PathName] = uigetfile('*.mat','Select the directory and file');
    cd(PathName);
    load(FileName, 'structvar');
    structvar=structvar;
    
    
    cd(currdir);


end