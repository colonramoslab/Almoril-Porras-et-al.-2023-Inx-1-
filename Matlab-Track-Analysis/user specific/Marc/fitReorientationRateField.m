function model = fitReorientationRateField(field_in_runs, field_in_reos, dt)
%function model = fitReorientationRateField(field_in_runs, field_in_reos, dt)
%function model = fitReorientationRateField(eset, fieldname)
%
%fits a model in which the reorientation rate is linearly proportional to
%field value (a + b*fieldvalue)

debug = true;

if (nargin < 2)
    disp ('fitReorientationRateField(field_in_runs, field_in_reos, dt) or  fitReorientationRateField(eset, fieldname)');
    model = [];
    return;
end

if (nargin == 2)
    eset = field_in_runs;
    fieldname = field_in_reos;
    field_in_runs = eset.gatherField(fieldname,'runs');
    field_in_reos = eset.gatherFromSubField('reorientation', fieldname, 'position', 'start');
    dt = mean(eset.gatherSubField('dr', 'interpTime'));
end

modelfun = @(params, val) params(1) + params(2)*val;

%minimizefun = @(x) -sum(log(modelfun(x,reoth))) + dt*sum(modelfun(x,allth)) + 1E12*(x(1) < abs(x(2)));
minimizefun = @(x) -logLikelihood(x, field_in_reos, field_in_runs, dt);

%----start fixing here ------%

fx = linspace(min(field_in_runs), max(field_in_runs), ceil(length(field_in_reos)/20));
ydat = histc(field_in_reos,fx)./histc(field_in_runs,fx)/dt;

x0 = polyfit(double(fx), double(ydat), 1);
x0 = x0([2 1]);
%x0 = [length(field_in_reos)/length(field_in_runs); 0];

%x0
%op = optimset('Algorithm','active-set');
%op.MaxFunEvals = 1000;
%op.GradObj = 'on';
%x = fmincon(minimizefun, x0, A, B, [], [], [], [],[],op);
%x
x = fminsearch(minimizefun, x0);
op.TolFun = 1.0000e-010;
op.TolX = 1.0000e-009;
x = fminsearch(minimizefun, x,op);
[ll,g,h] = logLikelihood(x, field_in_reos, field_in_runs, dt);

model.ll = ll;
model.grad = g;
model.hess = h;
model.params = x;

if (debug)
    plot (fx, ydat, 'b.-', fx,x(1)+x(2)*fx,'k--');
end

function [r,g,h] = rate(x, fv)
r = x(1) + x(2)*fv;
r(r < eps) = eps;
if (size(r,1) ~= 1)
    r = r';
end
if (nargout == 1)
    return;
end
g = zeros(2,length(fv));
g(1,:) = 0;
g(2,:) = fv;
h = zeros(2,3,length(fv));

function [ll, g, h] = logLikelihood(x, field_in_reos, field_in_runs, dt)
[rr,gr] = rate(x, field_in_reos);
[rn,gn] = rate(x, field_in_runs);
ll = sum(log(rr)) - dt * sum(rn);

if (nargout == 1)
    return;
end

g = sum(gr./[rr;rr],2) - dt*sum(gn,2);
rdiv = repmat(reshape(rr, 1, 1, []), [2 2]);

grtimesgr = zeros(2,2,length(gr));
grtimesgr(1,1,:) = gr(1,:).^2;
grtimesgr(2,2,:) = gr(2,:).^2;
grtimesgr(1,2,:) = gr(1,:).*gr(2,:);
grtimesgr(2,1,:) = grtimesgr(1,2,:);
whos grtimesgr
whos rdiv
h = sum(-grtimesgr./rdiv.^2,3);% + sum(hr./rdiv,3) - dt*sum(hn,3);
