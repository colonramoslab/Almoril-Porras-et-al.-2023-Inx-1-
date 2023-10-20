function [starts] = findFirstResponse(respCnt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
starts=NaN([size(respCnt,2),1]);
for ii=1:size(respCnt,2)
    foo=min(find(respCnt(:,ii)));
    if ~isempty(foo)
        starts(ii)=foo;
    end
end
end

