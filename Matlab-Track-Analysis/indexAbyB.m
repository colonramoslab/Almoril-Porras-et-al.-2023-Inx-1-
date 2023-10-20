function [outA,outB] = indexAbyB(inA,inB)
%UNTITLED6 Use values of B to create indices for values of A
% Inputs:
%   inA, a matrix with values for a particular feature of an object to be
%   ordered based on inB values
%   inB, a matrix with values for another feature that will be used to
%   create indices for outA
% Outputs:
%   outA,
%   outB, 

% issue: space is not descritized by pixel... Could multiply by 4sig figs
% or just round now...
inB=round(inB);
% linearize & remove NaN values for simplicity
inB=inB(~isnan(inB)); inA=inA(~isnan(inB));

% Initialize outB
% define index space, outB
x1=min(min(inB)); 
xend=max(max(inB));
outB=x1:xend;

% Intialize outA
% define depth of outA space
x=length(outB);
y=0;
for ii=outB
    ytest=sum(sum(inB==ii));
    y=max(y,ytest);
end
outA=nan([x,y]);

% Assign outA values from inA based on inB
for ii=1:x % length(outB)
    indVals=inA(inB==outB(ii));
    outA(ii,1:length(indVals))=indVals;
end   

end

