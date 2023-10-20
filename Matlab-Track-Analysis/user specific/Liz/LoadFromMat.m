function [eset] = LoadFromMat()
%function [] = LoadFromMat()
%
%Loads an experiment set from MAT files and includes utility to load from
%dialoge box

[fn, basedir] = uigetfile('*.mat','Select the .mat files you would like to use', 'MultiSelect', 'on');
if ~iscell(fn)
    fn = {fn};
end
basedirfn=[basedir fn{1}];
loadInds=strfind(basedirfn,'_');
loadFn=basedirfn(1:(loadInds(end-1)-1));
disp('loading...');
eset = ExperimentSet.fromMatFiles(loadFn);
end