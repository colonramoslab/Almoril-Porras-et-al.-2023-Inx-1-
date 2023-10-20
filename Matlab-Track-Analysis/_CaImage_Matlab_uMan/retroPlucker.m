function [ii] = retroPlucker(dirs)
%UNTITLED5 From input of old data (MetaMorph-derived), outputs variables
%for analysis pipeline developed for uManager-derived
%   

% Input should be cell array of pointers to old dataStructures

for ii=1:length(dirs)
    tDir=dirs{ii};
    [tsT,tsF, outDir] = retroPluckdStruct(tDir);

end

