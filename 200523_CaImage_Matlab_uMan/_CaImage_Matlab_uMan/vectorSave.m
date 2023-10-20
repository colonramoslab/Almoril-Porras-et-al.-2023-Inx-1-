function [ fStat, cnt ] = vectorSave( fignum, saveName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fStat=0;
cnt=0;
oWrite=1;

if ~ishandle(fignum)
    cnt='Need to input valid figure handle';
else
    set(fignum,'Units','Inches');
    pos=get(fignum,'Position');
    set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos(3),pos(4)]);
    fStat=1;
end

if isempty(saveName)
    saveName=fullfile(pwd,'figure.pdf');
end

[sdir,sn,ex]=fileparts(saveName);

if isempty(sdir)
    sdir=pwd;
    if isempty(ex)
        ex=('.pdf');
    end
    saveName=fullfile(sdir,strcat(sn,ex));
end

if ~oWrite
    while exist(saveName)
        cnt=cnt+1;
        warning('non-unique savename, modified to avoid overwrite');
        [sdir,sn,ex]=fileparts(saveName);
        saveName=fullfile(sdir,strcat(sn,num2str(cnt),ex));
    end
end
       
if fStat
    print(fignum,saveName,'-dpng','-r150','-painters');
end

end

