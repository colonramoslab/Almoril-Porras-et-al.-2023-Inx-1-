function [midline, gradmidline] = snakeToFindMidline(cpts, end1, end2, varargin)
%function midline = snakeToFindMidline(cpts, end1, end2, varargin)

%close all
nmidpts = 21;
alpha = 0.5;
beta = 0.1;

debug = true;
varagin = assignApplicable(varargin);


cpts = double(cpts);
end1 = double(end1);
end2 = double(end2);
%initialize midline guess
mid = splitOutline(cpts, end1, end2);

midlen = sum(sqrt(sum(diff(mid,[],2).^2)));

cpts = cpts/midlen;
mid = mid/midlen;

midold = mid;

mid = lowpass1D(mid, length(mid)/nmidpts);
mid(:,1) = midold(:,1);
mid(:,end) = midold(:,end);
mid = interp1(mid', linspace(1,length(mid), nmidpts))';
mid = resampleContour(mid, 'closed', false);

%create forcefield images
xaxis = linspace(1.1*min(cpts(1,:))-0.1*(max(cpts(1,:))), 1.1*max(cpts(1,:))-0.1*(min(cpts(1,:))), length(cpts));
yaxis = linspace(1.1*min(cpts(2,:))-0.1*(max(cpts(2,:))), 1.1*max(cpts(2,:))-0.1*(min(cpts(2,:))), length(cpts));


[xx,yy] = meshgrid(xaxis,yaxis);

in = inpolygon(xx(:), yy(:), cpts(1,:), cpts(2,:));

dt = DelaunayTri(cpts');
nn = dt.nearestNeighbor(xx(:), yy(:));
distim = sqrt((xx(:) - cpts(1,nn)').^2 + (yy(:) - cpts(2,nn)').^2);
distim(~in) = -10*distim(~in);
distim = reshape(distim, size(xx)) - max(distim(:));

[xforce, yforce] = imgradient(distim,5);
xstep = diff(xaxis(1:2));
ystep = diff(yaxis(1:2));
xforce = xforce/xstep;
yforce = yforce/ystep;

d1 = 0*sqrt((xx(:) - mid(1,1)).^2 + (yy(:) - mid(2,1)).^2) + sqrt((xx(:) - cpts(1,nn)').^2 + (yy(:) - cpts(2,nn)').^2);
d2 = 0*sqrt((xx(:) - mid(1,end)).^2 + (yy(:) - mid(2,end)).^2) + sqrt((xx(:) - cpts(1,nn)').^2 + (yy(:) - cpts(2,nn)').^2);
d1(~in) = 10*d1(~in);
d2(~in) = 10*d2(~in);

d1pot = 0.1*reshape(d1, size(distim)); %+ distim;
d2pot = 0.1*reshape(d2, size(distim)); %+ distim;
oldmid = mid;

%tdmat = diag(ones([1 nmidpts]) + alpha*tau, 0) - tau*alpha/2 * diag(ones([1 nmidpts-1]), 1)- tau*alpha/2 * diag(ones([1 nmidpts-1]), -1);
%whos tdmat
%figure(1);pcolor(tdmat); shading flat;

upot = []; ubend = []; utotal = []; uspr = [];
utotal = 1E9;

energyfun = @(x) energy(alpha, beta, [mid(:,1) x mid(:,end)], xaxis, yaxis, -distim);
energyfun2 = @(x) energy(alpha, beta, x, xaxis, yaxis, -distim);
gradenergyfn = @(x) gradOfEnergy(alpha, beta, x, xaxis, yaxis, -xforce, -yforce);

[utotal, uspr, ubend, upot] =  energy(alpha, beta, mid, xaxis, yaxis, -distim);
mid = conjugateGradient(energyfun2, gradenergyfn, mid, 'maxiter', 1000);

[utotalnew, usprnew, ubendnew, upotnew] =  energy(alpha, beta, mid, xaxis, yaxis, -distim);

clf(); pcolor (xaxis, yaxis, -distim); shading flat; hold on; plot(oldmid(1,:), oldmid(2,:), 'y.-', mid(1,:), mid(2,:), 'c.-'); hold off

disp(['du = ' num2str(utotalnew-utotal)]);
disp(['duspr = ' num2str(usprnew-uspr)]);
disp(['dubend = ' num2str(ubendnew-ubend)]);
disp(['dupot = ' num2str(upotnew - upot)]);
gx = gradenergyfn(mid);
midline = mid*midlen;
gradmidline = gx;
return


function dx = gradOfEnergy(alpha, beta, x, ux, uy, dux, duy)
n = length(x);
dcharge = alpha * gradchargeenergy(x, false);
dbend = beta * gradbendenergy(x,false);
fx = interp2(ux,uy,dux,x(1,:),x(2,:), '*linear', NaN);
fx(~isfinite(fx)) = mean(x(1,:)) - x(1,~isfinite(fx));
fy = interp2(ux,uy,duy,x(1,:),x(2,:), '*linear', NaN);
fy(~isfinite(fy)) = mean(x(2,:)) - x(2,~isfinite(fy));
dx = dcharge + dbend + [fx;fy]/n;


function [utotal, ucharge, ubend, upot] = energy(alpha, beta, x, ux, uy, u)
ucharge = alpha * chargeenergy(x, false);
ubend = beta * bendenergy(x, false);
upot = mean(interp2(ux,uy,u,x(1,:),x(2,:), '*linear', max(u(:)*100)));
utotal = ucharge + ubend + upot;% + upot1 + upot2;



