function [right, left] = skew_normal_pdfOfdThetaBreakdown (x,thetaIn, dTheta)
%function [right, left] = skew_normal_pdfOfdThetaBreakdown (x,thetaIn,
%dTheta)
% P(dt | theta) = (0.5 + g) * S(U, sigma, skew, dt) + (0.5 - g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];theta_0 = x(8);

theta_0 = x(7);
tt = thetaIn - theta_0;

g = x(3)*sin(tt);
u = x(1) - x(2)*cos(tt);
skew = x(4) - x(5)*cos(tt);
sigma = x(6);



%function pdf = skewNormalPDF(x, u, s, a)

right = (0.5 - g).*skewNormalPDF (dTheta, u, sigma, skew);
left =  (0.5 + g).*skewNormalPDF (-dTheta, u, sigma, skew);
end
