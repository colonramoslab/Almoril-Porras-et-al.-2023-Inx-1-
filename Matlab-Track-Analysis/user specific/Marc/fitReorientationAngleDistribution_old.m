function model = fitReorientationAngleDistribution(thetaIn, dTheta, theta_0)
%function model = FitReorientationAngleDistribution(thetaIn, dTheta, theta_0)
%
%assume P(dTheta) = (0.5 + g(theta)) N(u + f(theta), sigma) + (0.5 -
%g(theta)) N (u - f(theta), sigma)
%
%g(theta) = A + B(sin(theta - theta_0)) 
%f(theta) = C - D(cos(theta - theta_0))
%
%parameters to fit are B, C, D, sigma, theta_0:  set A, u to 0 (no
%asymmetry)
if (~exist('theta_0', 'var'))
    minimizefn = @(x) -logLikelihood(thetaIn, dTheta, 0, x(1), x(2), x(3), 0, x(4), x(5)) + 1E15*(abs(x(1)) > 0.5 | abs(x(3)) > abs(x(2)));
    theta_0 = atan2(mean(sin(thetaIn)), mean(cos(thetaIn)));
    %x0 = [0 0.2 pi/4 pi/8 0 pi/3 theta_0]; 
    x0 = [0.2 pi/4 pi/8 pi/3 theta_0]'; 
else
    minimizefn = @(x) -logLikelihood(thetaIn, dTheta, 0, x(1), x(2), x(3), 0, x(4), theta_0) + 1E15*(abs(x(1)) > 0.5 | abs(x(3)) > abs(x(2)));
%    theta_0 = atan2(mean(sin(thetaIn)), mean(cos(thetaIn)));
    %x0 = [0 0.2 pi/4 pi/8 0 pi/3 theta_0]; 
    x0 = [0.2 pi/4 pi/8 pi/3]'; 
end
for j = 1:2
    clear op;
    x = fminsearch(minimizefn, x0);
    op.TolFun = 1.0000e-010;
    op.TolX = 1.0000e-009;
    if (x(2) < 0)
        x(3) = -x(3);
        x(2) = -x(2);
        x(1) = -x(1);
    end
    x = fminsearch(minimizefn, x, op);
    if (x(3) < 0)
        x(5) = mod(x(3)+pi, 2*pi);
        x(3) = -x(3);
        x(1) = -x(1);
    end
    x = fminsearch(minimizefn, x, op);

    op.Display = 'off';
    op.DiffMinChange = 1E-24;
    op.FinDiffType = 'central';
    op.LargeScale = 'off';
    
    [x,fval1,exitflag,output,grad,hess1] = fminunc(minimizefn, x, op);
    hess = hessian(minimizefn, x);
    hess1
    hess
    if (all(isreal(hess)) && all(isreal(x)) && all(isfinite(x(:))) && all(isfinite(hess(:))))
        break;
    end
    inds = 1:length(x0);
    x0 = [0 pi/2 0 std(dTheta) 0]';
    x0 = x0(inds);
end
%{
for qq = 1:10
    try
        [v,d] = eig(hess);
        if (any(d < 0))
            disp ('not at local minimum');
            break
        end
        if ~all(isreal(hess))
            hess = hessian(minimizefn, x);
        end
        s = zeros(size(d)); s(logical(eye(size(s)))) = 0.1./sqrt(d(logical(eye(size(s)))));
        s = s.*randn(size(s));
        dx = sum(v*s,2);
        [xnew,fval,exitflag,output,grad,hess] = fminunc(minimizefn, x + dx, op);
        if (all(isfinite(xnew)) && all(isreal(xnew)) && all(isreal(hess(:))) && minimizefn(x) < minimizefn(xnew))
            x = xnew;
        end
    catch me
        continue;
    end
end
%}
%x
%fval
model.hessian = -hessian(minimizefn, x);
model.grad = -gradest(minimizefn, x);
model.params = x;
%model.pdfOfdTheta = @(x, thetaIn, dTheta) pdfOfdTheta(thetaIn, dTheta, x(1), x(2), x(3), x(4), x(5), x(6), x(7));
model.pdfOfdTheta = @(x, thetaIn, dTheta) pdfOfdTheta(thetaIn, dTheta, 0, x(1), x(2), x(3), 0, x(4), x(5));
model.paramkey = 'ratioAmplitude, offsetMean, offsetAmplitude, gaussSigma, theta_0';
model.ratioAmplitude = x(1);
model.offsetMean = x(2);
model.offsetAmplitude = x(3);
model.gaussSigma = x(4);
if (length(x) >= 5)
    model.theta_0 = mod(x(5) + pi, 2*pi) - pi;
else
    model.theta_0 = theta_0;
end
    
model.directionIndex = 2*model.ratioAmplitude.*model.offsetMean / sqrt(model.offsetMean.^2 + model.gaussSigma.^2);
model.magnitudeIndex = model.offsetAmplitude./sqrt(model.offsetMean.^2 + model.gaussSigma.^2);
model.reorientationMeansVsTheta = @(x,thetaIn) reorientationMeansVsTheta (thetaIn, x(1), x(2), x(3), x(4), x(5));
end

function p = pdfOfdTheta (thetaIn, dTheta, A, B, C, D, u, sigma, theta_0)
g = A + B*sin(thetaIn-theta_0);
f = C - D*cos(thetaIn-theta_0);

p = (0.5 + g).*normpdf(dTheta, u - f, sigma) + (0.5 - g).*normpdf(dTheta, u + f, sigma);
end

function [meandt, meanabsdt] = reorientationMeansVsTheta (thetaIn,  B, C, D, sigma, theta_0)
g = B*sin(thetaIn-theta_0);
u = C - D*cos(thetaIn-theta_0);
meandt = -2.*g.*u;
meanabsdt = sigma.*exp(-u.^2./(2*sigma.^2)) * sqrt(2/pi) + u.*erf(u./(sigma.*sqrt(2)));
end

function ll = logLikelihood (thetaIn, dTheta, A, B, C, D, u, sigma, theta_0)


ll = sum(log(pdfOfdTheta (thetaIn, dTheta, A, B, C, D, u, sigma, theta_0)));


end