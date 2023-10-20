function [xt] = xTzero(fh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% figure(fh);
xt=get(gca,'xticklabels');
xti=str2num(xt{1});
for ii=1:length(xt)
    xt{ii}=num2str(str2num(xt{ii})-xti);
end
set(gca,'xticklabels',xt)
end

