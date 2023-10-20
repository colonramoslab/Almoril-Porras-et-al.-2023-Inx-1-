

% kruskall-wallis anova, AFD & AIY threshold
groups={'AFD15','AFD20','AFD25','','','AIY15','AIY20','AIY25'};
g={}; x=[];
for i=1:3
    t=allTset(i,:);
    t=t(~isnan(t));
    x=[x,t];
    gt=repmat(groups(i),[1,numel(t)]);
    g=[g, gt]; 
end
[p,tbl,stats] = kruskalwallis(x,g);



%% non-parametric ranksum test
q1=qAIY15;
q2=qAIY20;
x=q1.Amp(~isnan(q1.Amp));
y=q2.Amp(~isnan(q2.Amp));
[p,h,stats] = ranksum(x,y); 


%% Fisher's exact test for dichotomous variable 
x = table([3;1],[6;7],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'})
[h,p,stats] = fishertest(x)