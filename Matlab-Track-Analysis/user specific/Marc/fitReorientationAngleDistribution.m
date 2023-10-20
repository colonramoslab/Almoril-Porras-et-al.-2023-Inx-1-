function model = fitReorientationAngleDistribution(thetaIn, dTheta, fixedValues, startValues, doTimeConsumingCalcs)
%function model = fitReorientationAngleDistribution(thetaIn, dTheta,
%fixedValues, startValues)
%
% P(dt | theta) = (0.5 - g) * S(U, sigma, skew, dt) + (0.5 + g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];
%to hold a quantity fixed, pass fixedValues with field set
%e.g. fixedValues.sigma = 0.3;



%P(dt | theta) = a*N(0,sigma1,dt) + (1/2 - a/2 + g) * S(U, sigma2, skew, dt) +
%    (1/2 - a/2 - g) * S(U, sigma2, skew, -dt)
%
%a(theta) = A + B*cos(theta - theta_0); (|B| < A);
%g(theta) = C(sin(theta - theta_0)) 
%other fit parameters, U, sigma1, sigma2, skew, theta_0 
%
%model.params = [A, B, C, U, sigma1, sigma2, skew, theta_0];
%to hold a quantity fixed, pass fixedValues with field set
%e.g. fixedValues.sigma1 = 0.3;
debug = false;
existsAndDefault('fixedValues', []);
existsAndDefault('startValues', []);
existsAndDefault('doTimeConsumingCalcs', true);
parameterNames = {'A', 'B', 'C', 'D', 'E', 'sigma', 'theta_0'};

A0 = 0.1; %0 <= A [x1]
B0 = 0; %|B| < A [x2]
C0 = 0.1; %|C| < 0.5 [x3]
D0 = 1; % 0 <= D; [x4]
E0 = 0; %|E| <= D; [x5]
sigma_0 = std(dTheta); % 0 <= sigma [x6]
theta_0_0 = atan2(mean(sin(thetaIn)), mean(cos(thetaIn))); %[x7]

%constraints:
%-1 * x1 < 0;
%1*x2 -1*x1 < 0;
%-1*x2 -1*x1 < 0;
%1 * x3 < 0.5;
%-1 * x3 < 0.5;
%-1*x4 < 0;
%1*x5 - x4 < 0;
%-1*x5 - x4 < 0;
%-1*x6 < 0;
B = [0 -eps -eps 0.5 0.5 0 0 0 0]';
x0 = [A0, B0, C0, D0, E0, sigma_0, theta_0_0]';

A = zeros(length(B), length(x0));
A(1,1) = -1;
A(2,1:2) = [-1 1];
A(3,1:2) = [-1 -1];
A(4,3) = -1;
A(5,3) = 1;
A(6,4) = -1;
A(7,4:5) = [-1 1];
A(8,4:5) = [-1 -1];
A(9,6) = -1;

minimizefn = @(x) -logLikelihood(x, thetaIn, dTheta);

Aeq = [];
Beq = [];

if (~isempty(startValues) && isstruct(startValues))
    fn = fieldnames(startValues);
    fnvalid = intersect(fn, parameterNames);
    fninvalid = setdiff(fn, parameterNames);
    if (~isempty(fninvalid))
        disp ('the following start parameters have invalid names:');
        disp (fninvalid);
        disp ('valid start parameter names are:');
        disp (parameterNames);
    end
    
    for j = 1:length(fnvalid)
        ind = strcmp(fnvalid{j}, parameterNames);
        x0(ind) = startValues.(fnvalid{j});
    end
end
if (~isempty(fixedValues) && isstruct(fixedValues))
    fn = fieldnames(fixedValues);
    fnvalid = intersect(fn, parameterNames);
    fninvalid = setdiff(fn, parameterNames);
    if (~isempty(fninvalid))
        disp ('the following fixed parameters have invalid names:');
        disp (fninvalid);
        disp ('valid fixed parameter names are:');
        disp (parameterNames);
    end
    
    Aeq = zeros(length(fnvalid), length(x0));
    Beq = zeros(length(fnvalid),1);
    for j = 1:length(fnvalid)
        ind = strcmp(fnvalid{j}, parameterNames);
        Aeq(j,ind) = 1;
        Beq(j) = fixedValues.(fnvalid{j});
        x0(ind) = fixedValues.(fnvalid{j});
    end
    
end

x = x0;
op.TolFun = 1.0000e-018;
op.TolX = 1.0000e-15;
op.Display = 'off';
op.DiffMinChange = 1E-24;
op.DiffMaxChange =  0.05;%000
op.FinDiffType = 'central';
op.LargeScale = 'off';
op.Algorithm = 'active-set';
op.MaxFunEvals = 6400;

if (any(A*x > B))
    disp (['equations violated: ']); disp( num2str(find(A*x > B)));
    disp ('initial condition violates constraints');
end


[x,FVAL,EXITFLAG,OUTPUT,LAMBDA,grad,hess] = fmincon(minimizefn, x, A, B, Aeq, Beq, [], [], [], op);%(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (debug)
    OUTPUT.message
    plotFit(x, thetaIn, dTheta);
end

model.logLikelihood = -FVAL;
model.minfun = minimizefn;
%{
model.hessian = -hess;
model.grad = -grad;
%}
if (doTimeConsumingCalcs)
    model.hessian = -hessian(minimizefn, x);
    model.grad = -gradest(minimizefn, x);
else
    model.hessian = [];
    model.grad = [];
end
model.params = x;
model.pdfOfdTheta = @skew_normal_pdfOfdTheta;% @(x, thetaIn, dTheta) pdfOfdTheta(x,thetaIn, dTheta);
model.fixedValues = fixedValues;
model.pdfOfdThetaBreakdown = @skew_normal_pdfOfdThetaBreakdown; %@(x, thetaIn, dTheta) pdfOfdThetaBreakdown(x,thetaIn, dTheta);
%model.plotFit = @(x, thetaIn, dTheta) plotFit(x, thetaIn, dTheta);
model.paramkey = parameterNames;
model.reorientationMeansVsTheta = @skew_normal_model_reorientationMeansVsTheta;

model.param_confidenceLevel = 0.95;
params_ci = zeros(2, length(x));
if (doTimeConsumingCalcs)
    for j = 1:length(x)
        [minv maxv] = confidenceRange(x, model.hessian, @(val) val(j,:), model.param_confidenceLevel);
        params_ci(:,j) = [minv maxv];
    end
end
model.params_ci = params_ci;

%{
model.ratioAmplitude = x(3);
model.offsetMean = x(4);
model.offsetAmplitude = 0;
model.gaussSigma = x(5);
if (length(x) >= 6)
    model.theta_0 = mod(x(6) + pi, 2*pi) - pi;
else
    model.theta_0 = theta_0;
end
    
model.directionIndex = 2*model.ratioAmplitude.*model.offsetMean / sqrt(model.offsetMean.^2 + model.gaussSigma.^2);
model.magnitudeIndex = model.offsetAmplitude./sqrt(model.offsetMean.^2 + model.gaussSigma.^2);
, x(1), x(2), x(3), x(4), 0, x(5), x(6));
%}
end

function [right, left] = pdfOfdThetaBreakdown (x,thetaIn, dTheta)
% P(dt | theta) = (0.5 + g) * S(U, sigma, skew, dt) + (0.5 - g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];theta_0 = x(8);

[right,left] = skew_normal_pdfOfdThetaBreakdown(x, thetaIn, dTheta);
% 
% theta_0 = x(7);
% tt = thetaIn - theta_0;
% 
% g = x(3)*sin(tt);
% u = x(1) - x(2)*cos(tt);
% skew = x(4) - x(5)*cos(tt);
% sigma = x(6);
% 
% 
% 
% %function pdf = skewNormalPDF(x, u, s, a)
% 
% right = (0.5 - g).*skewNormalPDF (dTheta, u, sigma, skew);
% left =  (0.5 + g).*skewNormalPDF (-dTheta, u, sigma, skew);
 end

function p = pdfOfdTheta (x,thetaIn, dTheta)
    p = skew_normal_pdfOfdTheta(x, thetaIn, dTheta);
end
% 
% [r,l] =  pdfOfdThetaBreakdown (x,thetaIn, dTheta);
% p = r + l;
% [r2,l2] = pdfOfdThetaBreakdown(x, thetaIn, 2*pi+dTheta);
% p = p + r2 + l2;
% [r2,l2] = pdfOfdThetaBreakdown(x, thetaIn, -2*pi+dTheta);
% p = p + r2 + l2;
% 
% if (any (p < 0))
%     %disp (['negative probability from x = ' num2str(x','%10.5e\t')]);
%   %  pause
% end
% if (any (~isfinite(p)))
%    % disp (['non-finite probability from x = ' num2str(x','%10.5e\t')]);
%    % pause
% end
% p(~isfinite(p)) = eps;
% p(p < 0) = eps;
% end


function [ll, g, h] = logLikelihood  (x,thetaIn, dTheta)
if (any(~isfinite(x)))
    ll = -Inf;
    return;
end

ll = sum(log(pdfOfdTheta  (x,thetaIn, dTheta)));

if (nargout > 1)
    g = gradest(@(x)  sum(log(pdfOfdTheta  (x,thetaIn, dTheta))), x);
end
if (nargout > 2)
    h = hessian(@(x)  sum(log(pdfOfdTheta  (x,thetaIn, dTheta))), x);
end
    

end

function stop = plotFit(x, thetaIn, dTheta)
    
    dtx = deg2rad(-170:20:170);
    dtxf = (-180:1:180);
    theta = [0 -pi pi/2 -pi/2];
    for j = 1:4
        inds = cos(theta(j)-thetaIn) > 1/sqrt(2);
        subplot(2,2,j); bar (rad2deg(dtx), hist(dTheta(inds), dtx)/nnz(inds)/diff(dtx(1:2)));
        [r, l] = pdfOfdThetaBreakdown(x, theta(j),deg2rad(dtxf));
        hold on; plot (dtxf, r+l, 'k-', dtxf, l, 'k--', dtxf,r, 'k--', 'LineWidth', 2); hold off
    end
%    disp(['x = ' num2str(x','%10.5e\t')]);
 %   disp(['ll = ' num2str(logLikelihood(x, thetaIn, dTheta))]);
    stop = false;
end