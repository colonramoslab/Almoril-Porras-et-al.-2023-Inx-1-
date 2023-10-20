function newpt = rethreshold (mra, pt, varargin)
%newpt = rethreshold (mra, pt, varargin)
%rethresholds the current point to produce
%a new MaggotTrackPoint
%updates the threshold contour targetArea area loc and cov 

targetArea = mra.targetArea;
scale = mra.contourScale;
debug = mra.debug;
varargin = assignApplicable(varargin);

newpt = MaggotTrackPoint (pt);
if (scale > 1)
    im2 = imresize(blurim(double(pt.imData), 1), scale, 'bicubic');
else
    scale = 1;
    im2 = double(pt.imData);
end
[bwim, thresh] = thresholdToTargetArea(im2, targetArea*scale^2);



bwim = centralRegion(bwim);
bwim = imfill(bwim, 'holes');

newpt.threshold = thresh;
newpt.targetArea = targetArea;
newpt.area = bwarea(bwim)/(scale^2);

%im = zeros(size(pt.imData));
%im(bwim) = pt.imData(bwim);

p = regionprops(bwim, im2, 'WeightedCentroid');

x = ((p.WeightedCentroid(1,:)) - [1 1])/scale + [1 1];
newpt.loc = (x + double(pt.imOffset) - [1 1])';

ctr = findContour (bwim);
ctr = ctr(:,1:(end-1));

scaledim = double(im2)./thresh;
s = [1 2 1; 0 0 0; -1 -2 -1];
[xd,yd] = imgradient(scaledim, 1);
energyim =  scale.^2*(-xd.^2-yd.^2);
ctr = marcsnake(1, 1, [], [], energyim, [3 1], ctr);

ctr = (ctr - ones(size(ctr)))/scale + ones(size(ctr));

offset = double(reshape(pt.imOffset,[2 1])) - [1;1];
newpt.contour = ctr + repmat(offset, 1, length(ctr));
try 
    cc = mra.track.expt.camcalinfo;
    if(~isempty(cc))
        newpt.contour = cc.realPtsFromCamPts(newpt.contour);
        newpt.loc = cc.realPtsFromCamPts(newpt.loc);
    end
catch
end
    

if debug
    figure(1);
    newpt.drawTrackImage();
end


function bwim = centralRegion (bwim)
%finds only the region that occupies 
%the center of the image
bwl = bwlabel(bwim);
sz = size(bwl);
sz = round(sz/2);
num = bwl(sz(1),sz(2));
if (num ~= 0)
    bwim = (bwl == num);
else
    p = regionprops(bwl, 'area');
    [a, I] = max([p.Area]);
    bwim = bwl == I;
end

function bwim = onlyLargest (bwim)
%function bwim = onlyLargest (bwim)
%finds only the largest region in bwim
%eliminates all other regions

bwl = bwlabel (bwim);
p = regionprops(bwl, 'area');

[a, I] = max([p.Area]);
bwim = bwl == I;

function contour = findContour (bwim)
%function contour = findContour (bwim)
%
%uses bwboundaries to find the boundary
%assumes you have one isolated blob with no holes

b = bwboundaries(bwim);
blah = b{1}';
contour = [blah(2,1:(end-1));blah(1,1:(end-1))];
