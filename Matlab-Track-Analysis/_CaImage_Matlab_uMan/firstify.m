function [outCnts] = firstify(respCnts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
outCnts={};

for ii=1:length(respCnts)
    rC=respCnts{ii};
    oC=zeros(size(rC));
    for jj=1:size(rC,2);
        fR=min(find(rC(:,jj)));
        if ~isempty(fR)
            oC(fR,jj)=1;
        end
    end
    outCnts{1,ii}=logical(oC);
end
end

