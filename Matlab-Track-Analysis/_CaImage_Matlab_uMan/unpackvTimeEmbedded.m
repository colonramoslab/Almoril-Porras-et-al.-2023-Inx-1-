function [vTimeSec,elapsedT,iDay] = unpackvTimeEmbedded(tR2)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
iDay=[];
vTime=cell([length(tR2),3]);
elapsedT=nan([length(tR2),1]);
% test=tR2{1,1};
% k=strfind(test,'"Time": "');
% iDay=test(k+9:k+18);
% % patch b/c of issue with embedded metadata reading, " are differently
% % encoded...
% if isempty(k)
%     k=strfind(test,'Time');
%     iDay=test(k(4)+7:k(4)+16);
% end

for ii=1:length(tR2)
    temp=tR2{ii};
%     tR3=strsplit(tR2{ii},',');
    % goal is elapsed time (sometimes 85)
    eP= strfind(temp,'ElapsedTime');
    elapsedT(ii)=str2double(temp(eP+16:eP+19));
    tP = strfind(temp,'"Time":');
    timeInfo= temp(tP+19:tP+26);
    tR7=(strsplit(timeInfo,':'));
    vTime(ii,:)=tR7;
    if isempty(iDay)
        iDay=temp(tP+8:tP+17);
    end
end

vTime=vTime(find(~cellfun(@isempty,vTime)));
vTime=str2double(vTime);
vTime=reshape(vTime,[length(vTime)/3,3]);

% Convert to common sec-based frame
vTime(:,1)=vTime(:,1)*60*60;
vTime(:,2)=vTime(:,2)*60;
vTimeSec=sum(vTime,2);
elapsedT=[0;elapsedT];
vTimeSec=[vTimeSec(1)-1;vTimeSec];



end

