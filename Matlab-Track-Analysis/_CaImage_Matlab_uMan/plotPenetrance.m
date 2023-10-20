function [ fh ] = plotPenetrance(dMat, varargin)
%UNTITLED11 Summary of this function goes here
%   dMat, columns are groups; rows are samples

[~, gCnt]=size(dMat);
groupLabels={};
compSets = []; % set of all pair-wise statistical comparisons, e.g. [1,2;1,3] for reference against two groups
mcc=1; % correct with bonferroni
alphThresh=0.05;

varargin=assignApplicable(varargin);

if isempty(groupLabels)
    for ii=1:gCnt
        groupLabels{ii}=num2str(ii);
    end
end

%% plot with binomial error,
fh=figure();
hold on;

totalG=sum(~isnan(dMat));
resp=sum(and(~isnan(dMat),dMat>0));
fResp=resp./totalG;
fError=sqrt((fResp.*(1-fResp))./totalG);
h1=errorbar(fResp,fError,'lineStyle','none');
h2=bar(fResp,'FaceColor',rand([1,3]));

if ~isempty(compSets)
    if mcc
        alphThresh=alphThresh/size(compSets,1);
    end
    
    for ii=1:size(compSets,1)
        
        % Fisher's exact test for dichotomous variable
        samp=compSets(ii,:);
        rowNames=groupLabels;
        varNames={'respond','silent'};
        x = table(resp(samp)',(totalG(samp)-resp(samp))','VariableNames',varNames, 'rowNames',rowNames(samp));
        [~,p,~] = fishertest(x);
        
        if p<alphThresh
            yL=get(gca,'ylim');
            line([compSets(ii,1),compSets(ii,2)],[yL(2),yL(2)]);
            text(mean([compSets(ii,1),compSets(ii,2)]),1.05*yL(2), '*','fontsize',25,'fontName','arial');
            set(gca,'ylim',[yL(1),1.1*yL(2)])
        end
    end
end


x=1:length(groupLabels);


ylabel('Fraction per group');
set(gca,'xtick',x,'xticklabels',groupLabels);
set(gca,'fontName','arial','fontSize',25);

end

