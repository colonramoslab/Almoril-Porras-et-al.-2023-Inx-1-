function [dstim, zoomcorners] = drawZoomBox(srcim, dstim, subimcent, startsize, newcent, endsize)
%function [dstim, zoomcorners] = drawZoomBox(srcim, dstim, subimcent, newcent, startsize, endsize)
%
%takes a subimage from srcim centered at subimcent (x,y) with dimensions
%(width x height), scales it to be endsize, then draws it into dstim 
%centered at newcent
%
%zoomcorners returns a cell {x,y} with the points on the corners of the new
%zoom rectangle, s.t. plot(zoomcorners{:}) will draw a nice box around the
%destination box

yindssrc = round((subimcent(2) - startsize(2)/2 - 0.5) + (1:startsize(2)));
xindssrc = round((subimcent(1) - startsize(1)/2 - 0.5) + (1:startsize(1)));

xindssrc = xindssrc(xindssrc > 0 & xindssrc <= size(srcim, 2));
yindssrc = yindssrc(yindssrc > 0 & yindssrc <= size(srcim, 1));


yindsdst = round((newcent(2) - endsize(2)/2 - 0.5) + (1:endsize(2)));
xindsdst = round((newcent(1) - endsize(1)/2 - 0.5) + (1:endsize(1)));


xindsdst = xindsdst(xindsdst > 0 & xindsdst <= size(dstim, 2));
yindsdst = yindsdst(yindsdst > 0 & yindsdst <= size(dstim, 1));

endsize = [length(yindsdst) length(xindsdst)];

subim = srcim(yindssrc,xindssrc);

dstim(yindsdst,xindsdst) = imresize(subim, endsize);

x1 = min(xindsdst);
x2 = max(xindsdst);
y1 = min(yindsdst);
y2 = max(yindsdst);

zoomcorners = {[x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1]};