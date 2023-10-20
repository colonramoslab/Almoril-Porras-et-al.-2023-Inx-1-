function [fh] = matchAxes(fh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
h=get(fh,'children');

ah2=h(2);
ah3=h(3);
pI2=ah2.Position;
pI3=ah3.Position;
pI3(1)=pI2(1); pI3(3)=pI2(3);
set(ah2,'Position',pI3);

end

