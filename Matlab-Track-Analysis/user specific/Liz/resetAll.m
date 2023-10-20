% resets all of the variables that check if you have completed an analysis
% stage in TemporalLightAnalysis so that you can load another experiment
% form scratch using the same program even if you have one open already on
% the workspace

names=fieldnames(ao);
 for j=1:length(names)
     ao.(names{j})=false;
     ao.resetAll=true;
 end