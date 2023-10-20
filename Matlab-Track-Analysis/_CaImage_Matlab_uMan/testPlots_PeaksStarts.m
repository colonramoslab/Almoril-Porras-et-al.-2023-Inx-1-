% overlay peak

figure(); imagesc(fluor')
hold on
pMat=sPeaks;
for ii=1:size(pMat,2)
    tt=pMat(:,ii);
    tp=find(tt);
    for jj=1:length(tp)
        plot(tp(jj),ii,'linestyle','none','marker','o','markersize',5,'markerfacecolor','w');
%         for jj=1:length(stRep)
%             [~] = line('XData',stRep(jj),'YData',ii,...
%                 'Color','w','Marker','o', 'MarkerSize',10,'lineStyle','none')
    end
end

peakPoints=peakMat>0;
pMat=peakPoints;
for ii=1:size(pMat,2)
    tt=pMat(:,ii);
    tp=find(tt);
    for jj=1:length(tp)
        plot(tp(jj),ii,'linestyle','none','marker','o','markersize',5,'markerfacecolor','g');
%         for jj=1:length(stRep)
%             [~] = line('XData',stRep(jj),'YData',ii,...
%                 'Color','w','Marker','o', 'MarkerSize',10,'lineStyle','none')
    end
end