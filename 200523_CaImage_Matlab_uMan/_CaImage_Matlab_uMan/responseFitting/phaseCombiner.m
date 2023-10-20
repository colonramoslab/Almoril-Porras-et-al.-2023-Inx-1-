function [comboMat,comboTime] = phaseCombiner(p1Mat,p1Time,p2Mat,p2Time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~(size(p1Time)==size(p2Time))
    warning('Matrices are two different sizes');
end


% % Make time indexable integer. Doesn't work b/c x1 & x2 are held as sig
% figs with really small non-obvious decimal values... silly to fix
% x1= nanmin(nanmin(diff(p1Time))); % Smallest increment of time
% x2= nanmin(nanmin(diff(p2Time))); % Smallest increment of time
% 
% if ~(x1==x2)
%     warning('time on different scales for two matrices');
% end
% 
% x1=0.0001
% n=0;
% while (floor(x1*10^n)~=x1*10^n)
%     n=n+1;
% end

% Assume at 10hz, then make integers
p1TimeIndex=p1Time*10;
p1TimeIndex=int16(p1TimeIndex);
p2TimeIndex=p2Time*10;
p2TimeIndex=int16(p2TimeIndex);
rel1Values=~(p1TimeIndex==0);
rel2Values=~(p2TimeIndex==0);

% used these indices to fill these matrices:
comboMat=nan(size(p1Time));
comboTime=nan(size(p1Time));    % Unnecessary, but use as test.

for ii=1:size(comboTime,2)

    % Assign P1 to column ii position
    colInd=rel1Values(:,ii); % indices to use in this column
    relP1now=p1TimeIndex(colInd,ii); % values to assign in this column
    comboTime(relP1now,ii)=p1Time(colInd,ii); % Assign times (validation)
    comboMat(relP1now,ii)=p1Mat(colInd,ii); % Assign values
    
    % Assign P2 to column ii position
    colInd=rel2Values(:,ii); % indices to use in this column
    relP2now=p2TimeIndex(colInd,ii); % values to assign in this column
    comboTime(relP2now,ii)=p2Time(colInd,ii); % Assign times (validation)
    comboMat(relP2now,ii)=p2Mat(colInd,ii); % Assign values

    
end
    
end

