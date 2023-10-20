% PULL set



%% Load tsF, fluorescence, and tsT, temperature traces with conditions
%genos={wt, wt, che2KO, che2KO, AFDr, AFDr} 

% WT
dirs={};
dirs{1}=fullfile(pwd,'demoData\181009_DCR7450_1');
dirs{2}=fullfile(pwd,'demoData\181009_DCR7450_2');
[temp_Wt, fluor_Wt] = loadTS(dirs);

 
% eat-4 che-2 KO
dirs={};
dirs{1}=fullfile(pwd,'demoData\181009_che2pFLP_DCR7452');
dirs{2}=fullfile(pwd,'demoData\181009_che2pflp_DCR7525_1');
[temp_che2KO, fluor_che2KO] = loadTS(dirs);
 
 
% eat-4 che-2 KO with AFD rescue
dirs={};
dirs{1}=fullfile(pwd,'demoData\181009_AFDeat4_DCR7521');
dirs{2}=fullfile(pwd,'demoData\181009_AFDeat4_DCR7522_1');
[temp_AFDresc, fluor_AFDresc] = loadTS(dirs);