function [RR] = rrDeRespCnt(respCnt, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

windO=[];
winVect=[];
sampS=size(respCnt,2);
protS=size(respCnt,1);

RR=zeros([sampS,1]);

varargin=assignApplicable(varargin);

if isempty(windO)
    % entire protocol
    windO=[1:protS];
end

if isempty(winVect)
    winVect=zeros([protS,1]);
    winVect(windO)=1;
end

winVect=logical(winVect);

for ii=1:sampS
    rT=respCnt(:,ii);
    RR(ii)=(10*nansum(rT(winVect))./protS)'; %in Hz, sampling at 10Hz
end

end

