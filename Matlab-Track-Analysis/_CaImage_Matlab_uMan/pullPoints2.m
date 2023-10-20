function [fluorMat, tempMat, peakMat ] = pullPoints2( fluor, temp, sPeaks, sampWin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%figure(); hold on;


peakPoints=find(sPeaks);
pkCnt=0;

fluorMat=NaN([2*sampWin,length(peakPoints)]);
tempMat=NaN([2*sampWin,length(peakPoints)]);
peakMat=NaN(size(sPeaks));

for i=1:size(sPeaks,2)
    % going through each sample to ensure within sample constraints,
    % - after beginning of recording and before end of recording
    sampFluor=fluor(:,i);
    sampTemp=temp(:,1);
    sampPeaks=find(sPeaks(:,i));
    for j=1:length(sampPeaks)
        % which peak is being filled
        pkCnt=pkCnt+1;
        % which data points to include
        ind=(sampPeaks(j)-sampWin):(sampPeaks(j)+sampWin-1);
        % are these data points within sampling window?
        offBeg=find(ind>0,1,'first'); % min(find(ind>0));
        offEnd=find(ind<length(fluor),1,'last'); % max(find(ind<length(fluor)));
        indReal=ind(offBeg:offEnd);
        % assign values to output matrices
        fluorMat(offBeg:offEnd,pkCnt)=sampFluor(indReal);
        tempMat(offBeg:offEnd,pkCnt)=sampTemp(indReal);
        
        % Find peak value for each response within fluorMat window
        peakVal=nanmax(fluorMat(50:offEnd,pkCnt));
        peakPos=find(sampFluor==peakVal);
        % select only within sample, Only important if many equivalent
        % points within sample, especially if large
        if length(peakPos)>1
            ppTemp=[];
            for kk=1:length(peakPos)
                if isempty(ppTemp)
                    ppTemp=indReal(peakPos(kk)==indReal);
                end
            end
            peakPos=ppTemp;
        end
        
        
        % may need to select only within ind(offBeg:offEnd)
        peakMat(peakPos,i)=peakVal;
        % testPlots_PeaksStarts.m
        
    end
end

