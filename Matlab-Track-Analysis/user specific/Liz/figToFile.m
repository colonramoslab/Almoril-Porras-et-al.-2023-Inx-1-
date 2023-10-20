function [] = figToFile(fstub)
%
%function [] = figToFile(fstub)
%this function prints all open figures to the file fstub with the name
%"figure#"

%create save directory if it doens't exist already
if exist(fstub,'file')~=7
    mkdir(fstub)
end

%print all open figures to file
figHandles=findobj('Type','figure');
count=1;
for n=1:length(figHandles)
    figure(figHandles(n));
    print('-depsc',[fstub 'figure' num2str(count)]);
    count=count+1;
end

end