function [y,J] = cosineForFitting(x,xdata)
%function [y,J] = cosineForFitting(x,xdata)

x(2) = x(2)/length(xdata); 
y = x(1)*cos(x(2)*xdata + x(3));
J(:,1) = cos(x(2)*xdata + x(3));
J(:,2) = -x(1)*xdata/length(xdata).*sin(x(2)*xdata + x(3));
J(:,3) = -x(1)*sin(x(2)*xdata + x(3));
%J = -J;