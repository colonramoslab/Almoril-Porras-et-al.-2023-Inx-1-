function [vTimeSec,elapsedT,iDay] = imfExtract(imFull)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%NEED QC in this for metadata format

imData=imfinfo(imFull);

% "Time" buried in UnknownTags!!
tR2=cell([length(imData),1]);
for ii=2:length(imData)
    tR2{ii}=imData(ii).UnknownTags.Value(2:end-1);
end
tR2=tR2(2:end)';

[vTimeSec,elapsedT,iDay] = unpackvTimeEmbedded(tR2);
%
%     mD1=cell([length(imData),1]);
%     mD2=nan([length(imData),1]);
%     mD2(1)=0;
%
%     for ii=2:length(imData)
%         t=imData(ii).UnknownTags.Value;
%         t=strsplit(t,',');
%         t1=t{48};
%         t1=strsplit(t1,'"');
%         mD1{ii}=t1{4};
%         t2=t{85};
%         t2=strsplit(t2,':');
%         mD2(ii)=str2num(t2{2});
%     end


end

