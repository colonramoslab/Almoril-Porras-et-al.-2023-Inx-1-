% This script is located in:
% ..._CaImage_Matlab_uMan\responseFitting\exampleScript.m


%% To adapt for your analysis:
% 1. Change dirs1 & dirs2 cell pointers to locations of data for groups 1 & 2
% 2. Make sure x1 & x2 encompass the AFD firing window, which will change if
% you are using a different protocol from 
% 3. Change the legValues to explain relevant differences/conditions. NOTE.
% These are used to save the file by default
% 4. You can also change the save directory using saveDir.


%% Tc=15, cat-1

% parent directory where data are located
dDir1='C:\Users\jshha\Dropbox\ColonRamosLab\Manuscripts\Manuscript_2019_GluNP\NTscreen_CaImaging\Joon'; 

% save to current working directory, or specify here.
saveDir=pwd;

% Add location with new function to path 
funPath='C:\Users\jshha\Dropbox\ColonRamosLab\Matlab\CaImaging_uMan\_CaImage_Matlab_uMan';
addpath(genpath(funPath));
    % FYI, New functions here: 'C:\Users\jshha\Dropbox\ColonRamosLab\Matlab\CaImaging_uMan\_CaImage_Matlab_uMan\responseFitting';


% Values for legend in figures
legValues{1}='Wt-15-AFD';
legValues{2}='cat1-15-AFD';

% Analysis window for protocol
x1 =1555;
x2 =1815;
incFrames{1}=[x1:x2];
incFrames{2}=incFrames{1}; % use the same window in both conditions.

% % (olaIs17 background control Tc=15C) DCR3056 olaIs17
dirs1={};
dirs1{1}=fullfile(dDir1,'190601_NT_mutants\190601f_olaIs17_Tc15_1');
dirs1{2}=fullfile(dDir1,'190601_NT_mutants\190601l_olaIs17_Tc15_1');
dirs1{3}=fullfile(dDir1,'190605_NT_mutants\190605b_olaIs17_Tc15_1');

% % (cat-1 Tc=15C) DCR7436 olaIs17; cat-1(e1111)
dirs2={};
dirs2{1}=fullfile(dDir1,'190528_NT_mutants\190528_run5_cat-1_Tc15_1');
dirs2{2}=fullfile(dDir1,'190528_NT_mutants\190528_run6_cat-1_Tc15_1');


[comboMat1,comboTime1,comboMat2,comboTime2] = responseCompare(dirs1,dirs2,legValues, 'incFrames', incFrames, 'saveDir', saveDir); 