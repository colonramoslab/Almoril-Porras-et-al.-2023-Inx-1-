function [outputArg1,outputArg2] = plotMeanResponse(respMats,tempMats, varargin)
%UNTITLED11 Summary of this function goes here
%   dMat, columns are groups; rows are samples

[~, gCnt]=size(dMat);
groupLabels={};
compSets = []; % set of all pair-wise statistical comparisons, e.g. [1,2;1,3] for reference against two groups
alphThresh=0.05;

varargin=assignApplicable(varargin);


for ii=1:length(respMats)
    
    tM=tempMats{ii};
    rM=respMats{ii};
    
    % gross mean features
    figure(); hold on;
    yyaxis left; ll=plot(mean(tM,2));
    yyaxis right; rl=plot(mean(rM,2));
   
    % plot each resonse as function of features, temperature & dTemp
    
    tR= mean(tM(45:55,:),1); % mean of the middle 11 frames
    
    
    
end

