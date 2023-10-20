function [xc,yc,structvar] = detectCheckerBoard(xaxis, yaxis, im, varargin)
%function [xc,yc, structvar] = detectCheckerBoard(xaxis, yaxis, im, varargin)
%
%   finds the corners of a checkerboard image, with some help
%   xaxis, yaxis are the axes for the image (i.e. pcolor(axis,yaxis,im)
%   makes sense)
%   im must be a grayscale image (size nrowsxncolsx1)
%optional key/value paris
%"nrows" number of rows of the checkerboard (y-direction)
%"ncols" number of columns of the checkerboard (x-direction)
%
%xc,yc are coordinates of corners
%structvar follows the rules of structvar, which are
%
%about structvar
%the first rule of structvar is you do not talk about structvar
%the second rule of structvar is you do not talk about structvar
%the third rule of structvar is that it is a Nx2 array of points,
%specifying rectangles
%the fourth rule of structvar is the lower left (in image coordinates; upper left in xy coords)
%corner comes first, then the other 3 points are specified in
%counterclockwise (in image coords; clockwise in xy coords) order
%the fifth rule of structvar is you do not talk about structvar
%
%structvar specifies the boundaries of light squares

nrows = [];
ncols = [];
varargin = assignApplicable(varargin);
if (isempty(xaxis))
    xaxis = 1:size(im,2);
end
if (isempty(yaxis))
    yaxis = 1:size(im,1);
end



figure(10); clf(10); 
pcolor(xaxis, yaxis, im); shading flat; colormap gray

if isempty(nrows)
    nrows = input ('how many rows');
end

if isempty(ncols)
    ncols = input ('how many columns');
end


%look for corner locations by convolving with checker pattern

%construct kernel
ksize = [size(im,1)/(3*nrows),  size(im,2)/(3*ncols)];
%make ksize odd
ksize = 2 * round(ksize/2)+1;

kernel = zeros(ksize);
kernel(1:floor(ksize(1)/2),1:floor(ksize(2)/2)) = 1;
kernel(ceil(ksize(1)/2 + 1):end,1:floor(ksize(2)/2)) = -1;
kernel(ceil(ksize(1)/2 + 1):end,ceil(ksize(2)/2 + 1):end) = 1;
kernel(1:floor(ksize(1)/2),ceil(ksize(2)/2 + 1):end) = -1;

%convolve & discard outer region 
cim = zeros(size(im));
cimall = conv2(im, kernel, 'same');
yinds = round(ksize(1)/3):round(size(cim,1) - ksize(1)/3);
xinds = round(ksize(2)/3):round(size(cim,2) - ksize(2)/3);
cim(yinds,xinds) = cimall(yinds,xinds);
%pcolor(xaxis,yaxis,cim); shading flat; colormap jet;

%look for local extrema

%do a 2nd deriv filter
sd = zeros(size(cim));
tempim = conv2(gaussKernel(ksize(1)/18), gaussKernel(ksize(2)/18), cim, 'same') - conv2(gaussKernel(ksize(1)/6), gaussKernel(ksize(2)/6), cim, 'same');
yinds = round(ksize(1)/2):round(size(cim,1) - ksize(1)/2);
xinds = round(ksize(2)/2):round(size(cim,2) - ksize(2)/2);
sd(yinds, xinds) = tempim(yinds,xinds);
sdlg = abs(sd) > percentile(abs(sd), 0.99);

figure(11); clf(11);
pcolor(xaxis, yaxis, sd); shading flat; colormap jet;
figure(12); clf(12);
pcolor(xaxis, yaxis, sd.*sdlg); shading flat; colormap jet;


%find the locations of the corners
stats = regionprops(imdilate(sdlg,strel('disk',10, 0)), abs(cim), 'WeightedCentroid');
loc = [stats.WeightedCentroid];
[xl,yl] = localmaxima(abs(cim), loc(1:2:end), loc(2:2:end),max(ksize)/4);
%return back to coordinate system specified by xaxis,yaxis
%the extra shift of 1/2 was determined empirically
xl = interp1(xaxis, xl + 1/2, 'linear');
yl = interp1(yaxis, yl + 1/2, 'linear');

needsfixing = length(xl) ~= nrows*ncols;
figure(10); clf(10);
hold off
pcolor(xaxis,yaxis,im); shading flat; colormap gray; hold on; 
plot (xl,yl, 'yo', xl, yl, 'y+', 'MarkerSize', 20,'LineWidth',2);

if (~needsfixing)
        s = input('Do you want to add/delete corners y/[n]','s');
    needsfixing = strcmpi(s(1),'y');
else
    disp ('number of corners detected does not match nrows*ncols; you must fix');
end

if (needsfixing)
   
    while(1)
        figure(10); clf(10);
        hold off
        pcolor(xaxis,yaxis,im); shading flat; colormap gray; hold on; 
        plot (xl,yl, 'yo', 'MarkerSize', 10,'LineWidth',2);
        title ({['n points = ' num2str(length(xl)) ' ; target = ' num2str(nrows*ncols)], 'please click on corners to remove, then hit enter'});
        [x,y] = getpts;
        reminds = zeros([1 length(x)]);
        for j = 1:length(x)
            ds = (xl-x(j)).^2 + (yl-y(j)).^2;
            [blah,I] = min(ds);
            reminds(j) = I;
        end
        try 
            plot (xl(reminds), yl(reminds), 'rx', 'MarkerSize', 20,'LineWidth',2);
            s = input('accept subtractions? y/[n]','s');
            if (strcmpi(s(1),'y'))
                keep = true(size(xl));
                keep(reminds) = false;
                xl = xl(keep);
                yl = yl(keep); 
                if (length(xl) <= nrows*ncols)
                    break;
                end
            end
        catch me
            reminds
            
            disp('had a problem, try again, error follows:');
            disp(me.getReport);
        end
    end
    while(length(xl) < nrows*ncols)
        figure(10); clf(10);
        hold off
        pcolor(xaxis,yaxis,im); shading flat; colormap gray; hold on; 
        plot (xl,yl, 'yo', 'MarkerSize', 10,'LineWidth',2);
        title ({['n points = ' num2str(length(xl)) ' ; target = ' num2str(nrows*ncols)], 'please click on markers to add, then hit enter'});
        [x,y] = getpts;
        
        if (~isempty(x))
            x = interp1(xaxis, 1:length(xaxis), x);
            y = interp1(yaxis, 1:length(yaxis), y);
            [x,y] = localmaxima(abs(cimall), x, y);
            x = interp1(xaxis, x + 1/2, 'linear');
            y = interp1(yaxis, y + 1/2, 'linear');

            
            figure(10); hold on
            plot (x, y, 'go', 'MarkerSize', 20,'LineWidth',2);
            s = input('accept additions? y/[n]','s');
            if (strcmpi(s(1),'y'))
                xl = [xl x']; %#ok<AGROW>
                yl = [yl y'];  %#ok<AGROW>
                break;
            end
        end
    end
end
xc = xl;
yc = yl;

%turn corners into rectangles
[rs1,rs2] = corners2rectangles(xc,yc,ncols);

%find the value at the center of the first 2 rectangles
pt1 = round(mean(rs1(1:4,:), 1));
pt2 = round(mean(rs2(1:4,:), 2));
val1 = interp2(xaxis,yaxis,im, pt1(1), pt1(2));
val2 = interp2(xaxis,yaxis,im, pt2(1), pt2(2));

%the brighter rectangle gets turned into structvar
if (val1 > val2)
    structvar = rs1;
else
    structvar = rs2;
end




function [x,y] = localmaxima(im, x, y, nsize, maxiters, tol)
%function [x,y] = localmaxima(im, x, y, maxiters, tol)
%uses steepest ascent to move to local maxima;  
%this is really only a good idea if you're already close

debug = false;

existsAndDefault('nsize', 10);
existsAndDefault('maxiters', 100);
existsAndDefault('tol', 0.001);

xd = conv2(gaussKernel(nsize), dgausskernel(nsize), im, 'same');
yd = conv2(dgausskernel(nsize), gaussKernel(nsize), im, 'same');
ds = sqrt(xd.^2 + yd.^2);
xd = xd./ds;
yd = yd./ds;

xx = 1:size(im,2);
yy = 1:size(im,1);

xnew = x;
ynew = y;
nv = interp2(xx,yy,im,xnew,ynew);
if (debug)
    figure(1); clf(1); imagesc(im); 
    hold on
end
for j = 1:maxiters
    ov = nv;
    xdd = interp2(xx,yy,xd,x,y,'*linear');
    ydd = interp2(xx,yy,yd,x,y,'*linear');
    for k = 1:length(x)
        xline = x(k) + xdd(k)*(0:0.01:1)*nsize;
        yline = y(k) + ydd(k)*(0:0.01:1)*nsize;
        valsonline = interp2(xx,yy,im,xline,yline,'*linear');
        [blah,I] = max(valsonline);
        xnew(k) = xline(I);
        ynew(k) = yline(I);
    end
    if (debug)
        for k = 1:length(x)
            plot ([x(k) xnew(k)], [y(k) ynew(k)], 'b-', 'LineWidth',2);
        end
    end
    x = xnew;
    y = ynew;
    nv = interp2(xx,yy,im,xnew,ynew);
    if all((nv - ov)./ov < tol)
        return
    end
end

    