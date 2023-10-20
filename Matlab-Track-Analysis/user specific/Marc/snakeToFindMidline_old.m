function midline = snakeToFindMidline(cpts, end1, end2, varargin)
%function midline = snakeToFindMidline(cpts, end1, end2, varargin)

%close all
nmidpts = 11;
alpha = 4;%-6E-5;
beta = 0.5;%40E-8;

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
%figure(1); subplot(1,2,1); pcolor(xaxis, yaxis, d1pot); shading flat; hold on;
%plot (mid(1,:), mid(2,:), 'k-', cpts(1,:), cpts(2,:), 'c-', 'LineWidth', 3); hold off
%subplot(1,2,2); pcolor(xaxis, yaxis, d2pot); shading flat; hold on
%plot (mid(1,:), mid(2,:), 'k-', cpts(1,:), cpts(2,:), 'c-','LineWidth', 3); hold off
%figure(1); subplot(1,2,1); pcolor(xaxis,yaxis, xforce); shading flat; colorbar vert; axis equal
%subplot(1,2,2); pcolor(xaxis,yaxis, yforce); shading flat; colorbar vert; axis equal
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
%mid = fminsearch(energyfun2, mid);
mid = conjugateGradient(energyfun2, gradenergyfn, mid);
%{
xstart = mid(:,2:end-1);
xfin = fminsearch(energyfun, xstart);
mid(:,2:end-1) = xfin;
l = [0 cumsum(sum(diff(mid,[],2).^2))];
mid = interp1(l, mid', linspace(0, max(l), length(mid)))';
xstart = mid(:,2:end-1);
xfin = fminsearch(energyfun, xstart);
mid(:,2:end-1) = xfin;
l = [0 cumsum(sum(diff(mid,[],2).^2))];
mid = interp1(l, mid', linspace(0, max(l), length(mid)))';
%}

[utotalnew, usprnew, ubendnew, upotnew] =  energy(alpha, beta, mid, xaxis, yaxis, -distim);

disp(['du = ' num2str(utotalnew-utotal)]);
disp(['duspr = ' num2str(usprnew-uspr)]);
disp(['dubend = ' num2str(ubendnew-ubend)]);
disp(['dupot = ' num2str(upotnew - upot)]);

midline = mid*midlen;
return


function dx = gradOfEnergy(alpha, beta, x, ux, uy, dux, duy)

n = length(x);
dx = diff(x,[],2);
dl = sqrt(sum(dx.^2));

dx1 = [[0;0] -dx./([dl.^3;dl.^3])];
dx2 = [dx./([dl.^3;dl.^3]) [0;0]];
dspr = alpha * (dx1 + dx2)/(n^2);

dbend = n.^4*beta*[[0 0;0 0] diff(x,4,2) [0 0; 0 0]]; %use the 4th derivative and call it quits

fx = interp2(ux,uy,dux,x(1,:),x(2,:), '*linear', NaN);
fx(~isfinite(fx)) = mean(x(1,:)) - x(1,~isfinite(fx));
fy = interp2(ux,uy,duy,x(1,:),x(2,:), '*linear', NaN);
fy(~isfinite(fy)) = mean(x(2,:)) - x(2,~isfinite(fy));
dx = dspr + dbend + [fx;fy];


function [utotal, uspr, ubend, upot] = energy(alpha, beta, x, ux, uy, u)

n = length(x);
uspr = alpha*sum(1./sqrt(sum(diff(x,[],2).^2)))/(n^3);

%l = [0 cumsum(sum(diff(x,[],2).^2))];
%n = 100;
%x = interp1(l, x', linspace(0, max(l), n))';
dx = diff(x,[],2);
dls = sum(dx.^2);



that = dx./sqrt([dls;dls]);
cs = sum(diff(that, [],2).^2); %curvature squared
ubend = n*beta*sum(cs);
%ddx = deriv(dx, sigma);
%ubend = beta/n * sum((dx(1,:).*ddx(2,:) - dx(2,:).*dx(1,:)).^2./(dls.^3));
%ubend = beta*n^3*sum(sum(ddx.^2));
upot = sum(interp2(ux,uy,u,x(1,:),x(2,:), '*linear', max(u(:)*100)))/n;
%upot1 = interp2(ux,uy,uend1, x(1,1), x(2,1), '*linear',max(uend1(:)));
%upot2 = interp2(ux,uy,uend2, x(1,end), x(2,end), '*linear',max(uend2(:)));
utotal = uspr + ubend + upot;% + upot1 + upot2;



