function [newpt,gradmidline]  = bestMidline(track, pt, varargin)
%function midline = snakeToFindMidline(cpts, end1, end2, varargin)

%close all
nmidpts = 21;
alpha = 4;
beta = 1.4;
debug = true;
distpower = 2;
distoffset = 0.2;
varagin = assignApplicable(varargin);

if (~isa (pt, 'TrackPoint'))
    pt = track.pt(pt);
end
cc = track.expt.camcalinfo;

cpts = double(pt.contour);
sp = double(pt.spine);
if (~isempty(cc))
    cpts = cc.camPtsFromRealPts(cpts);
    sp = cc.camPtsFromRealPts(sp);
end
cpts = cpts - repmat(double(pt.imOffset') - [1;1], [1 length(cpts)]);
sp = sp - repmat(double(pt.imOffset') - [1;1], [1 length(sp)]);

im = (double(pt.imData)./double(pt.threshold));
[xd, yd] = imgradient(im, 1);
edgeim = sqrt(xd.^2 + yd.^2);
im = blurim(im, 2);


mid = sp;

midlen = sum(sqrt(sum(diff(mid,[],2).^2)));

[ends,ccurv] = findPointyEnds(cpts);
[~, c1, c2] = splitOutline(cpts, ends(1), ends(2));
cwidth = median(sqrt(sum((c1-c2).^2)));
%ccurv = 1:length(cpts);
curvim = zeros(size(im)) + max(abs(ccurv));
inds = sub2ind(size(im), round(cpts(2,:)), round(cpts(1,:)));
curvim(inds) = -ccurv*midlen;
curvim = blurim(imerode(curvim, ones(3)),0.5);

d1 = sum((mid(:, [1 1]) - cpts(:, ends)).^2);
if (d1(1) < d1(2))
    mid(:,1) = cpts(:,ends(1));
    mid(:,end) = cpts(:, ends(2));
else
    mid(:,1) = cpts(:,ends(2));
    mid(:,end) = cpts(:, ends(1));
end


cpts = cpts/midlen;

mid = mid/midlen;
xaxis = (1:size(im,2))/midlen;
yaxis = (1:size(im,1))/midlen;
midold = mid;

mid = lowpass1D(mid, length(mid)/nmidpts);
mid(:,1) = midold(:,1);
mid(:,end) = midold(:,end);
mid = interp1(mid', linspace(1,length(mid), nmidpts))';
mid = resampleContour(mid, 'closed', false);


[xx,yy] = meshgrid(xaxis,yaxis);

in = inpolygon(xx(:), yy(:), cpts(1,:), cpts(2,:));

dt = DelaunayTri(cpts');
nn = dt.nearestNeighbor(xx(:), yy(:));

distim = sqrt((xx(:) - cpts(1,nn)').^2 + (yy(:) - cpts(2,nn)').^2)*midlen/cwidth;%.^(distpower/2);
%endenergyim = reshape(100*distim.^2, size(xx));%-20*edgeim;
distim(~in) = -distim(~in);

distim = -(distim - max(distim(in)))./distoffset;
distim = reshape(distim, size(xx));
distim = distim.^(distpower);

%{
figure(2); clf(2)
pcolor(xaxis, yaxis, -distim); shading interp; caxis([0 mdi]); colorbar vert
hold on; plot (cpts(1,:), cpts(2,:), 'k-', 'LineWidth', 2); hold off
pause
%}
firstenergyim = -1*(distim < max(1, (.2/distoffset)^distpower));


energyim = blurim(distim,1);% -blurim(distim,1);
endenergyim = blurim(curvim,1);

xstep = diff(xaxis(1:2));
ystep = diff(yaxis(1:2));

[xgrad, ygrad] = imgradient(energyim,0.5);
xgrad = xgrad/xstep;
ygrad = ygrad/ystep;

[endxgrad, endygrad] = imgradient(endenergyim,0.5);
endxgrad = endxgrad/xstep;
endygrad = endygrad/ystep;

[firstxgrad, firstygrad] = imgradient(firstenergyim,0.5);
firstxgrad = firstxgrad/xstep;
firstygrad = firstygrad/ystep;


oldmid = mid;
alpha1 = 100;
beta1 = 1;

firstenergyfun = @(x) firstEnergy(alpha1, beta1, x, xaxis, yaxis, firstenergyim);
gradfirstenergyfun = @(x) gradfirstEnergyFixedEnds (alpha1, beta1, x, xaxis, yaxis, firstxgrad, firstygrad);
[utotal, uspr, ubend, upot] = firstenergyfun(mid);
mid = conjugateGradient(firstenergyfun, gradfirstenergyfun, mid, 'maxiter', 100, 'ftol', 1e-6); %pull it straight
[utotalnew, usprnew, ubendnew, upotnew] =  firstenergyfun(mid);
if(true && debug)
    disp(['du = ' num2str(utotalnew-utotal)]);
    disp(['duspr = ' num2str(usprnew-uspr)]);
    disp(['dubend = ' num2str(ubendnew-ubend)]);
    disp(['dupot = ' num2str(upotnew - upot)]);
    figure(3); clf(3);
    pcolor(xaxis, yaxis, firstenergyim); shading flat; colormap gray; axis equal; colorbar vert; hold on;
    plot (oldmid(1,:), oldmid(2,:), 'y.-', mid(1,:), mid(2,:), 'c.-', 'LineWidth', 2, 'MarkerSize', 20); 
    [gx, gspring, gbend, gpot] = gradfirstenergyfun(mid);
    gx = -gx/10;
    gspring = -gspring/10;
    gbend = -gbend/10;
    gpot = -gpot/10;
    %plot (mid(1,:), mid(2,:), 'b.-', 'LineWidth', 2, 'MarkerSize', 20); hold on;        
    quiver (mid(1,:), mid(2,:), gx(1,:), gx(2,:), 0, 'c');
    quiver (mid(1,:), mid(2,:), gspring(1,:), gspring(2,:), 0, 'm');
    quiver (mid(1,:), mid(2,:), gbend(1,:), gbend(2,:), 0, 'y');
    quiver (mid(1,:), mid(2,:), gpot(1,:), gpot(2,:), 0, 'r');
    hold off
end
oldmid = mid;
energyfun = @(x) energy(alpha, beta, x, xaxis, yaxis, energyim, endenergyim);
energyfun1 = @(x) energy(alpha, beta, x, xaxis, yaxis, energyim, zeros(size(endenergyim)));
gradenergyfn1 = @(x) gradOfEnergy(alpha, beta, x, xaxis, yaxis, xgrad, ygrad, endxgrad, endygrad, true);
gradenergyfn2 = @(x) gradOfEnergy(alpha, beta, x, xaxis, yaxis, xgrad, ygrad, endxgrad, endygrad, false);

[utotal, uspr, ubend, upot] = energyfun(mid);% energy(alpha, beta, mid, xaxis, yaxis, energyim);


mid = conjugateGradient(energyfun1, gradenergyfn1, mid, 'maxiter', 100, 'debug', false, 'conjugate', true, 'ftol', 1E-7);
%mid = conjugateGradient(energyfun, gradenergyfn2, mid, 'maxiter', 1000, 'debug', debug, 'conjugate', true, 'ftol', 1E-7);
mid = resampleContour(mid, 'closed', false);

if(debug)
    figure(1); clf(1);
    pcolor(xaxis, yaxis, energyim); shading interp; colormap jet; axis equal; colorbar vert;
    caxis([min(energyim(in)) max(energyim(in))]);
%    figure(2); clf(2);
 %   pcolor(xaxis, yaxis, endenergyim); shading flat; colormap jet; axis equal; colorbar vert;
    for mm = 1:1
        figure(mm); hold on;
        [gx, gcharge, gbend, gpot] = gradenergyfn1(mid);
        gx = -gx/10;
        gcharge = -gcharge/10;
        gbend = -gbend/10;
        gpot = -gpot/10;
        plot (oldmid(1,:), oldmid(2,:), 'y.-', mid(1,:), mid(2,:), 'b.-', 'LineWidth', 2, 'MarkerSize', 20); hold on;        
        quiver (mid(1,:), mid(2,:), gx(1,:), gx(2,:), 0, 'c');
        quiver (mid(1,:), mid(2,:), gcharge(1,:), gcharge(2,:), 0, 'm');
        quiver (mid(1,:), mid(2,:), gbend(1,:), gbend(2,:), 0, 'y');
        quiver (mid(1,:), mid(2,:), gpot(1,:), gpot(2,:), 0, 'r');
        [gx, gcharge, gbend, gpot] = gradenergyfn1(oldmid);
        gx = -gx/10;
        gcharge = -gcharge/10;
        gbend = -gbend/10;
        gpot = -gpot/10;
        quiver (oldmid(1,:), oldmid(2,:), gx(1,:), gx(2,:), 0, 'c');
        quiver (oldmid(1,:), oldmid(2,:), gcharge(1,:), gcharge(2,:), 0, 'm');
        quiver (oldmid(1,:), oldmid(2,:), gbend(1,:), gbend(2,:), 0, 'y');
        quiver (oldmid(1,:), oldmid(2,:), gpot(1,:), gpot(2,:), 0, 'r');
        hold off;
        legend ('im','oldctr','ctr', 'gx', 'charge', 'bend', 'potential');
        caxis([min(energyim(in)) max(energyim(in))]);
    end
    
end

[utotalnew, usprnew, ubendnew, upotnew] =  energyfun(mid); %energy(alpha, beta, mid, xaxis, yaxis, energyim);
gx = gradenergyfn1(mid);
if (false)
    clf(); pcolor (xaxis, yaxis, energyim); colormap jet; shading flat; hold on; plot(oldmid(1,:), oldmid(2,:), 'y.-', mid(1,:), mid(2,:), 'c.-'); quiver( mid(1,:), mid(2,:), gx(1,:), gx(2,:)); hold off
    colorbar vert; axis equal
end
if (debug)
    disp(['du = ' num2str(utotalnew-utotal)]);
    disp(['duspr = ' num2str(usprnew-uspr)]);
    disp(['dubend = ' num2str(ubendnew-ubend)]);
    disp(['dupot = ' num2str(upotnew - upot)]);
end
gx = gradenergyfn1(mid);
%mid = resampleContour(mid, false);
midline = mid*midlen;
gradmidline = gx;
whos midline

newpt = pt;
if (~isempty(cc))
    midline = cc.realPtsFromCamPts(midline + repmat(double(pt.imOffset') - [1;1], [1, length(midline)]));
end
newpt.spine = midline;

return


function [dx, dcharge, dbend, dpot] = gradOfEnergy(alpha, beta, x, ux, uy, dux, duy, edux, eduy, fixends)
n = length(x);
dcharge = -alpha * gradconstforceenergy(x, false);
dbend = beta * gradbendenergy(x,false);
fx = interp2(ux,uy,dux,x(1,:),x(2,:), '*linear', NaN);
fx(1) = interp2(ux,uy,edux,x(1,1),x(2,1), '*linear', NaN);
fx(end) = interp2(ux,uy,edux,x(1,end),x(2,end), '*linear', NaN);
if (any(~isfinite(fx)))
    disp ('contour out of range!');
end
fx(~isfinite(fx)) = mean(x(1,:)) - x(1,~isfinite(fx));
fy = interp2(ux,uy,duy,x(1,:),x(2,:), '*linear', NaN);
fy(1) = interp2(ux,uy,eduy,x(1,1),x(2,1), '*linear', NaN);
fy(end) = interp2(ux,uy,eduy,x(1,end),x(2,end), '*linear', NaN);
fy(~isfinite(fy)) = mean(x(2,:)) - x(2,~isfinite(fy));
dpot = [fx;fy]/n;
dx = dcharge + dbend + dpot;
if (fixends)
    dx(:,1) = 0;
    dx(:,end) = 0;
    dcharge(:,1) = 0;
    dcharge(:,end) = 0;
    dbend(:,1) = 0;
    dbend(:,end) = 0;
    dpot(:,1) = 0;
    dpot(:,end) = 0;
end


dxmag = sqrt(sum(sum(dx.^2)));
dcharge = dcharge/dxmag;
dbend = dbend/dxmag;
dpot = dpot/dxmag;
dx = dx/dxmag;

function [utotal, ucharge, ubend, upot] = energy(alpha, beta, x, ux, uy, u, eu)
n = length(x);
ucharge = -alpha * constforceenergy(x, false);
ubend = beta * bendenergy(x, false);
upot = 1/n*(sum(interp2(ux,uy,u,x(1,2:end-1),x(2,2:end-1), '*linear', max(u(:)*100))) + sum(interp2(ux,uy,eu,x(1,[1 end]),x(2,[1 end]), '*linear', max(u(:)*100))));
utotal = ucharge + ubend + upot;% + upot1 + upot2;

function [utotal, uspring, ubend, upot] = firstEnergy(alpha, beta, x, ux, uy, u)
uspring = alpha*springEnergy(x, false);
ubend = beta * bendenergy(x, false);
upot = mean(interp2(ux,uy,u,x(1,:),x(2,:), '*linear', max(u(:)*100)));
utotal = uspring + upot + ubend;

function [dx, dspring, dbend, dpot] = gradfirstEnergyFixedEnds(alpha, beta, x, ux, uy, dux, duy)
n = length(x);
dspring = alpha * gradspringenergy(x, false);
dspring(:,1) = 0;
dspring(:,end) = 0;
dbend = beta * gradbendenergy(x,false);
dbend(:,1) = 0;
dbend(:,end) = 0;
fx = interp2(ux,uy,dux,x(1,:),x(2,:), '*linear', NaN);
fx(~isfinite(fx)) = mean(x(1,:)) - x(1,~isfinite(fx));
fy = interp2(ux,uy,duy,x(1,:),x(2,:), '*linear', NaN);
fy(~isfinite(fy)) = mean(x(2,:)) - x(2,~isfinite(fy));
dpot = [fx;fy]/n;
dpot(:,1) = 0;
dpot(:,end) = 0;
dx = dpot + dspring;

