function [vTimeSec,elapsedT,iDay] = mmExtract(imFull)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% display(imFull);
H=tiff_read_header(imFull);

fileTime = ones(size(H,2),5);

for i=1:size(H,2)
    sliceTime = H{1,i}.DateTime;
    %trim irrelevant details
    temp = strsplit(sliceTime, {' ', ':', '.'} );
    tempnum = ones(1,size(temp,2));
    for j=1:size(temp,2)
        tempnum(j)=str2double(temp{j});
    end
    fileTime(i,:)=tempnum;
end

iDay=num2str(fileTime(1,1)); % 2018-10-09
iDay=strcat(iDay(1:4),'-',iDay(5:6),'-',iDay(7:8));

vTime=fileTime(:,2:end);
vTime(:,1)=vTime(:,1)*60*60;
vTime(:,2)=vTime(:,2)*60;
vTime(:,4)=vTime(:,4)/1000;
vTimeSec=sum(vTime,2);
elapsedT=vTimeSec-vTimeSec(1);
end

