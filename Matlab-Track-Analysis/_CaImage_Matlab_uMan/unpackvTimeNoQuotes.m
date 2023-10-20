function [vTimeSec,elapsedT,iDay] = unpackvTimeEmbedded(tR2)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
vTime=cell([length(tR2),3]);
elapsedT=nan([length(tR2),1]);
test=tR2{1,1};
k=strfind(test,'"Time": "');
iDay=test(k+9:k+18);
% patch b/c of issue with embedded metadata reading, " are differently
% encoded...
if isempty(k)
    k=strfind(test,'Time');
    iDay=test(k(4)+7:k(4)+16);
end

for ii=1:length(tR2)
    tR3=strsplit(tR2{ii},',');
    eTemp=strsplit(tR3{85},':');
    elapsedT(ii)=str2double(eTemp(2));
    tR4=strsplit(tR3{48},'"');
    tR5=tR4{1,4};
    tR6=tR5(12:19);
    tR7=transpose(strsplit(tR6,':'));
    vTime(ii,:)=tR7;
end
vTime=vTime(find(~cellfun(@isempty,vTime)));
vTime=str2double(vTime);
vTime=reshape(vTime,[length(vTime)/3,3]);

% Convert to common sec-based frame
vTime(:,1)=vTime(:,1)*60*60;
vTime(:,2)=vTime(:,2)*60;
vTimeSec=sum(vTime,2);
elapsedT=elapsedT(4:end);
end

