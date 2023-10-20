function [resps] = findResp(gArray)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sampDim=2; % dimension that identifies individual samples/replicates
resps=nan([nansum(nansum(gArray)),1]);
rI=1;

for ii=1:size(gArray,sampDim)
    rE=rI+nansum(gArray(:,ii))-1;
    resps(rI:rE)=find(gArray(:,ii)==1);
    rI=rE+1;
end

