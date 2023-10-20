function [ fh ] = barFractionPerGroup(dMat, varargin)
%UNTITLED11 Summary of this function goes here
%   dMat, columns are groups; rows are samples

[~, gCnt]=size(dMat);
groupLabels={};

varargin=assignApplicable(varargin);

if isempty(groupLabels)
    for i=1:gCnt
        groupLabels{i}=num2str(i);
    end
end

fh=figure();
hold on;

groups=unique(dMat);
groups=groups(~isnan(groups));
groups=sort(groups,'descend');

top=ones([1,size(dMat,2)]);
totalG=sum(~isnan(dMat));

% lay down bars
for i=1:length(groups)
    bar(top,'FaceColor',rand([1,3]));
    fG=sum(dMat==groups(i))./totalG;
    top=top-fG; 
end

% add labels
top=ones([1,size(dMat,2)]);
for i=1:length(groups)
    fG=sum(dMat==groups(i))./totalG;
    x=1:length(fG);
    y=(2*top-fG)./2;
    labelTxt=num2str(groups(i));   
    text(x,y,labelTxt);
    top=top-fG; 
end

ylabel('Fraction per group');
set(gca,'xtick',x,'xticklabels',groupLabels);


end

