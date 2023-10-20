function [fluorMat, tempMat ] = pullPoints( fluor, temp, peakPoints, sampWin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%figure(); hold on;


fluorMat=NaN([2*sampWin,length(peakPoints)]);
tempMat=NaN([2*sampWin,length(peakPoints)]);
for i=1:length(peakPoints)
    % define sampling window
    ind=(peakPoints(i)-sampWin):(peakPoints(i)+sampWin-1);
    % trim indices to fit constraints of data, >0 & <length of recording.
    offBeg=find(ind(ind>0),1,'first'); % min(find(ind>0));
    offEnd=find(ind(ind<length(fluor)),1,'last'); % max(find(ind<length(fluor)));
    % assign data to output matrix
    fluorMat(offBeg:offEnd,i)=fluor(ind(offBeg:offEnd));
    tempMat(offBeg:offEnd,i)=temp(ind(offBeg:offEnd));

end
end

