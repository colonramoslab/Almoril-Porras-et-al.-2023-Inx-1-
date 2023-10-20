function model = fitReorientationRateTheta(eset, allth, reoth)
%function model = fitReorientationRateTheta(eset, allth, reoth)

deltaLogLikelihood = log(100);

if (~existsAndDefault('allth', []))
    allth = eset.gatherField('theta','runs');
end
if (~existsAndDefault('reoth', []))
   reoth = eset.gatherSubField('reorientation', 'prevDir');
end

dt = eset.expt(1).track(1).dr.interpTime;

modelfun = @(params, theta) params(1) - params(2)*cos(theta - params(3));

%minimizefun = @(x) -sum(log(modelfun(x,reoth))) + dt*sum(modelfun(x,allth)) + 1E12*(x(1) < abs(x(2)));
minimizefun = @(x) -logLikelihood(x, reoth, allth, dt)+ 1E12*(x(1) < abs(x(2)));
tx = linspace(0, 2*pi, 10);
hh = histc(mod(allth, 2*pi), tx);
[~,I] = max(hh);

x0 = [length(reoth)/(length(allth)*dt); length(reoth)/(length(allth)*dt)/2; tx(I)];
A = [-1 1 0; 0 0 1; 0 0 -1];
B = [0; 6*pi; 0];
%op = optimset('Algorithm','active-set');
%op.MaxFunEvals = 1000;
%op.GradObj = 'on';
%x = fmincon(minimizefun, x0, A, B, [], [], [], [],[],op);
%x
x = fminsearch(minimizefun, x0);
op.TolFun = 1.0000e-010;
op.TolX = 1.0000e-009;
x = fminsearch(minimizefun, x,op);
[ll,g,h] = logLikelihood(x, reoth, allth, dt);

model.ll = ll;
model.grad = g;
model.hess = h;
model.params = x;

params = x;

model.fun = str2func(func2str(modelfun)); % try to strip out extraneous crap from handle when saving later
model.params = x;
model.index = x(2) / x(1); %amplitude / mean
model.direction = mod(x(3) + pi, 2*pi) - pi;


%{
ll = sum(log(modelfun(params, reoth))) - dt*sum(modelfun(params,allth));
parallelfun = @(y,theta) repmat(y(:,1), [1 length(theta)]) -  repmat(y(:,2),[1 length(theta)]).*(cos(repmat(theta, size(y,1), 1) -repmat(y(:,3),[1 length(theta)])));
model.llreal = ll;



dll = @(y) sum(log(parallelfun(y, reoth)),2) - dt*sum((parallelfun(y, allth)),2) -ll;

dll(params')
simplex1 = linspace(0.5*x(1), 1.5*x(1), 21);
simplex2 = linspace(0.5*x(2), 1.5*x(2), 21);
simplex3 = x(3) + linspace(-pi/4,pi/4, 21);

params = repmat(x', length(simplex1), 1);
params(:,1) = simplex1;
dloglik1 = dll(params);
figure(1)
plot (simplex1-x(1), dloglik1);

minx1 = simplex1(max(1,find(dloglik1 > -2*deltaLogLikelihood, 1, 'first')- 1));
maxx1 = simplex1(min(length(simplex1),find(dloglik1 > -2*deltaLogLikelihood, 1, 'last')+ 1));

params = repmat(x', length(simplex2), 1);
params(:,2) = simplex2;
dloglik2 = dll(params);
minx2 = simplex2(max(1,find(dloglik2 > -2*deltaLogLikelihood, 1, 'first')- 1));
maxx2 = simplex2(min(length(simplex2),find(dloglik2 > -2*deltaLogLikelihood, 1, 'last')+ 1));
figure(2)
plot (simplex2-x(2), dloglik2);

params = repmat(x', length(simplex2), 1);
params(:,2) = simplex2;
dloglik2 = dll(params);
minx2 = simplex2(max(1,find(dloglik2 > -2*deltaLogLikelihood, 1, 'first')- 1));
maxx2 = simplex2(min(length(simplex2),find(dloglik2 > -2*deltaLogLikelihood, 1, 'last')+ 1));

params = repmat(x', length(simplex3), 1);
params(:,3) = simplex3;
dloglik3 = dll(params);
figure(3)
plot (simplex3-x(3), dloglik3);
pause(0.1);
minx3 = simplex3(max(1,find(dloglik3 > -2*deltaLogLikelihood, 1, 'first')- 1));
maxx3 = simplex3(min(length(simplex3),find(dloglik3 > -2*deltaLogLikelihood, 1, 'last')+ 1));
if (isempty(minx1))
    minx1 = 0.99*x(1);
    maxx1 = 1.01*x(1);
end

x1 = linspace(minx1,maxx1, 21);
x2 = linspace(minx2,maxx2, 21);
x3 = linspace(minx3,maxx3, 21);

deltacube = zeros(length(x2), length(x1), length(x3));
[xx1,xx2] = meshgrid(x1, x2);
for j = 1:length(x3)
    params = [xx1(:),xx2(:),repmat(x3(j), length(xx1(:)),1)];
    deltacube(:,:,j) = reshape(dll(params), length(x2), length(x1));
end

model.px1 = x1;
model.px2 = x2;
model.px3 = x3;
model.deltacube = deltacube;


try
    dx1 = x1-x(1);
    dx2 = x2-x(2);
    figure(4);[c,h]=contour(dx1, dx2, exp(deltacube(:,:,11)), 1./([2 5 10 20 50 100])); clabel(c,h);
    hold on
    hh = model.hess(1:2,1:2);
    bob = ([dx1' dx2']*hh*[dx1;dx2]);
    [c,h]=contour(dx1, dx2, exp(bob), 1./([2 5 10 20 50 100])); clabel(c,h);
    hold off
catch me
    disp(me.getReport);
end
%}

function [r,g,h] = rate(x, theta)
r = x(1) - x(2)*cos(theta - x(3));
if (size(r,1) ~= 1)
    r = r';
end
if (nargout == 1)
    return;
end
g = zeros(3,length(theta));
g(1,:) = 1;
g(2,:) = -cos(theta-x(3));
g(3,:) = x(2)*sin(x(3)-theta);
h = zeros(3,3,length(theta));
h(2,3,:) = sin(x(3)-theta);
h(3,2,:) = h(2,3,:);
h(3,3,:) = x(2)*cos(x(3)-theta);

function [ll, g, h] = logLikelihood(x, theta_reo, theta_no_reo, dt)
if (nargout == 1)
    rr = rate(x, theta_reo);
    rn = rate(x, theta_no_reo);
    ll = sum(log(rr)) - dt * sum(rn);
    return;
end
[rr,gr,hr] = rate(x, theta_reo);
[rn,gn,hn] = rate(x, theta_no_reo);
ll = sum(log(rr)) - dt * sum(rn);

g = sum(gr./[rr;rr;rr],2) - dt*sum(gn,2);
rdiv = repmat(reshape(rr, 1, 1, []), [3 3]);

grtimesgr = zeros(3,3,length(gr));
grtimesgr(1,1,:) = gr(1,:).^2;
grtimesgr(2,2,:) = gr(2,:).^2;
grtimesgr(3,3,:) = gr(3,:).^2;
grtimesgr(1,2,:) = gr(1,:).*gr(2,:);
grtimesgr(1,3,:) = gr(1,:).*gr(3,:);
grtimesgr(2,1,:) = grtimesgr(1,2,:);
grtimesgr(3,1,:) = grtimesgr(1,3,:);
grtimesgr(2,3,:) = gr(2,:).*gr(3,:);
grtimesgr(3,2,:) = grtimesgr(2,3,:);



h = sum(-grtimesgr./rdiv.^2,3) + sum(hr./rdiv,3) - dt*sum(hn,3);
