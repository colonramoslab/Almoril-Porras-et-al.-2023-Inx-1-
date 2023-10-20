function [vTimeSec,elapsedT,iDay] = unpackvTime(tR2)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
vTime=cell([length(tR2),3]);
elapsedT=nan([length(tR2),1]);
test=tR2{1,1};
k=strfind(test,'"Time": "');
iDay=test(k+9:k+18);


for i=1:length(tR2)
    tR3=strsplit(tR2{i},',');
    % extract elapsed time
    tK=strfind(tR3,'lapsed');
    kI=find(~(cellfun(@isempty,tK)));
    
    if ~isempty(kI)
    eTemp=strsplit(tR3{kI},':');
    elapsedT(i)=str2double(eTemp(2));
    % extract clock time
    tK=strfind(tR3,'"Time": "');
    kI=find(~(cellfun(@isempty,tK)));
    tR4=strsplit(tR3{kI},'"');
    tR5=tR4{1,4};
    tR6=tR5(12:19);
    tR7=transpose(strsplit(tR6,':'));
    vTime(i,:)=tR7;
    end
end
vTime=vTime(find(~cellfun(@isempty,vTime)));
vTime=str2double(vTime);
vTime=reshape(vTime,[length(vTime)/3,3]);

% Convert to common sec-based frame
vTime(:,1)=vTime(:,1)*60*60;
vTime(:,2)=vTime(:,2)*60;
vTimeSec=sum(vTime,2);
elapsedT=elapsedT(~isnan(elapsedT));
% elapsedT=elapsedT(4:end);
end

