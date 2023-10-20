function model = fitReorientationAngleDistributionWithHSInfo(prevDir, nextDir, numHS, varargin)
%function model = fitReorientationAngleDistributionWithHSInfo(prevDir, nextDir, numHS, varargin)

debug = false;
doTimeConsumingCalcs = true;
varargin = assignApplicable(varargin);

thetaIn = prevDir;
dt = diff(unwrap([prevDir;nextDir]));
%{
fixedValues.A = 0;
fixedValues.B = 0;
fixedValues.sigma1 = std(dt(numHS == 0));
%}
snp = @(x, xd) max(skewNormalPDF(xd, x(1), x(2), x(3)), eps);
skewnormfun = @(x) -sum(log(snp(x,abs(dt(numHS > 0)))));

x = [mean(abs(dt(numHS > 0)))/2, std(abs(dt(numHS > 0))), 2];

x = fminsearch(skewnormfun, x);
dtx = 0:10:180;

if (debug)
    figure(1); clf
    bar (dtx, hist(abs(dt(numHS > 0)), deg2rad(dtx))/(sum(numHS > 0) * deg2rad(diff(dtx(1:2))))); hold on;
    plot (dtx, skewNormalPDF(deg2rad(dtx), x(1), x(2), x(3)), 'r-', 'LineWidth', 2);
    pause(0.1);
end
% P(dt | theta) = (0.5 - g) * S(U, sigma, skew, dt) + (0.5 + g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];



startValues.A = x(1);
startValues.sigma = x(2);
startValues.D = x(3);

model = fitReorientationAngleDistribution(prevDir(numHS > 0), dt(numHS > 0), [], startValues, doTimeConsumingCalcs);
%{
for m = 1:2
    fn = model.paramkey;
    for j = 1:length(fn)
        startValues.(fn{j}) = model.params(j);
    end
    clear fixedValues
    fn = setdiff(fn, {'A', 'B', 'C'});%{'U', 'sigma1', 'sigma2', 'skew', 'theta_0'};
    for j = 1:length(fn)
        fixedValues.(fn{j}) = startValues.(fn{j});
    end

    model = fitReorientationAngleDistribution(prevDir, dt, fixedValues, startValues);
    if (debug)
        plotFit(model, prevDir, dt, numHS)
        model.logLikelihood
        pause(0.1);
    end
    
    fn = model.paramkey;
    for j = 1:length(fn)
        startValues.(fn{j}) = model.params(j);
    end
    clear fixedValues
    fn = setdiff(fn, {'U', 'sigma2', 'skew'});%removed sigma1
    for j = 1:length(fn)
        fixedValues.(fn{j}) = startValues.(fn{j});
    end

    model = fitReorientationAngleDistribution(prevDir, dt, fixedValues, startValues);
    if (debug)
        plotFit(model, prevDir, dt, numHS)
        model.logLikelihood
        pause(0.1);
    end

    fn = model.paramkey;
    for j = 1:length(fn)
        startValues.(fn{j}) = model.params(j);
    end
    fn = setdiff(fn, {'theta_0'});%fn = {'A', 'B', 'C', 'U', 'sigma1', 'sigma2', 'skew'};
    for j = 1:length(fn)
        fixedValues.(fn{j}) = startValues.(fn{j});
    end
    model = fitReorientationAngleDistribution(prevDir, dt,fixedValues, startValues);
    if (debug)
        plotFit(model, prevDir, dt, numHS)
        model.logLikelihood
        pause(0.1);
    end
    
    
    fn = model.paramkey;
    for j = 1:length(fn)
        startValues.(fn{j}) = model.params(j);
    end
    
    clear fixedValues
    fn = {'sigma1'};
    for j = 1:length(fn)
        fixedValues.(fn{j}) = startValues.(fn{j});
    end

    model = fitReorientationAngleDistribution(prevDir, dt, fixedValues, startValues);
    if (debug)
        plotFit(model, prevDir, dt, numHS)
        model.logLikelihood
        pause(0.1);
    end
    
end
%}
%{
fun = @(x) sum(log(model.pdfOfdTheta(x, prevDir(numHS > 0), dt(numHS > 0))));
model.grad = gradest(fun, model.params);
model.hessian = hessian(fun, model.params);
%}
%model.plotFit = @(model, prevDir, dt, numHS) plotFit(model,prevDir, dt, numHS);
if (debug)
    try
        plotFit(model, prevDir, dt, numHS)
    catch me
        disp(me.getReport);
    end
end

end



function plotFit(model, thetaIn, dTheta, nhs)
    x = model.params;
    bs = 10;
    dtx = deg2rad((-180+bs/2):bs:(180-bs/2));
%    dtx = deg2rad(-170:20:170);
    dtxf = (-180:1:180);
    theta = [0 -pi pi/2 -pi/2];
    for j = 1:4
        inds = cos(theta(j)-thetaIn) > 1/sqrt(2);
        %h1 = hist(dTheta(inds & nhs == 0), dtx)/nnz(inds)/diff(dtx(1:2));
        h2 = hist(dTheta(inds & nhs > 0), dtx)/nnz(inds & nhs > 0)/diff(dtx(1:2));
        subplot(2,2,j); bar (rad2deg(dtx), h2);%[h1;h2]','stacked');
        [r, l] = model.pdfOfdThetaBreakdown(x, theta(j),deg2rad(dtxf));
        hold on; plot (dtxf, r+l, 'k-', dtxf, l, 'k--', dtxf,r, 'k--', 'LineWidth', 2); hold off
        title (['prev dir = ', num2str(rad2deg(theta(j)))]);
        sum((r+l)*deg2rad(diff(dtxf(1:2))))
        sum(h2*diff(dtx(1:2)))
    end
%    ttl = {['x = ' num2str(x','%10.5e\t')], ['ll = ' num2str(model.logLikelihood)]};
 %   title(ttl);
   
end
