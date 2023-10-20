function [thresh, area] = optimumThreshold(im,varargin)
%function [thresh, area] = optimumThreshold(im,varargin)
%
%optional args
%minThresh = 1;
%maxThresh = 255;
%minArea = 5;
%maxArea = 500;
%sigma = 1;
%showGraphs = false;

minThresh = 1;
maxThresh = 255;
minArea = 5;
maxArea = 500;
sigma = 1;
showGraphs = false;
assignApplicable(varargin);

m = double(ceil(max(im(:))));
if (m < maxThresh) 
    maxThresh = m;
end

m = double(ceil(min(im(:))));
if (m > minThresh)
    minThresh = m;
end

d = maxThresh - minThresh + 1;
area = zeros([1 d]);

for j = 1:d
    t = minThresh - 1 + j;
    area(j) = sum(im(:) >= t);
end
    
dv = deriv(log(area), sigma);
dv2 = deriv(dv, 1);
mint2 = find((area < maxArea) & (dv2 > 0), 1, 'first');
maxt2 = find((area > minArea) & (dv2 < 0), 1, 'last');
if (showGraphs)
    figure(1);
    semilogy(minThresh-1+(1:d),area,'b-');
    ylabel('log area'); xlabel('thresh');
    figure(2);
    plot (minThresh-1+(1:d),dv, minThresh-1+(mint2:maxt2),dv(mint2:maxt2),'r-');
    ylabel('D(log area)'); xlabel('thresh');
end
[y,I] = max(dv(mint2:(maxt2)));
area = area(I+mint2);
thresh = minThresh+mint2+I;