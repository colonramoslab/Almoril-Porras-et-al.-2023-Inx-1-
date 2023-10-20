function [outCell] = comboReady(gCell)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

outCell={};

% find sample dimension
% assumption, longest dimension.
sDim=find(size(gCell)==max(size(gCell)));

for ii=1:size(gCell,sDim)
    
    % find dimension of individual cell, should be longest
    dDim=find(size(gCell{ii})==max(size(gCell{ii})));
    if dDim==1
    outCell{ii}=gCell{ii}';
    else
        outCell{ii}==gCell{ii};
    end

end

