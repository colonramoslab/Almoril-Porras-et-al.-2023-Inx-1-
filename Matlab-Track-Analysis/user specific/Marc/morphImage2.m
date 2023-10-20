function imout = morphImage2(src, srcxaxis, srcyaxis, dstx, dsty, srcx, srcy, dstxaxis, dstyaxis)
%function imout = morphImage2(src, srcxaxis, srcyaxis, dstx, dsty, srcx,
%srcy, dstxaxis, dstyaxis)
%
%image is a picture taken in camera coordinates
%camerax,cameray <---> dstx,dsty are corresponding camera points and dst
%points
%dstxaxis,dstyaxis are axes in dst coordinates that the image is to be
%morphed to;  note that this may be any subset of the locations spanned by
%dstx,dsty

if (size(dstx,1) == 1)
    dstx = dstx';
end
if (size(dsty,1) == 1)
    dsty = dsty';
end
if (size(srcx,1) == 1)
    srcx = srcx';
end
if (size(srcy,1) == 1)
    srcy = srcy';
end


if (isempty(srcxaxis))
    srcxaxis = 1:size(src,2);
end
if (isempty(srcyaxis))
    srcyaxis = 1:size(srcy,1);
end

if (isempty(dstxaxis))
    dstxaxis = linspace(min(dstx), max(dstx), length(srcxaxis));
end
if (isempty(dstyaxis))
    dstyaxis = linspace(min(dsty), max(dsty), length(srcyaxis));
end


[rxp,ryp] = meshgrid(dstxaxis,dstyaxis);


[srcx, srcy, dstx, dsty] = guessOutsideHull (srcx, srcy, dstx, dsty, dstxaxis, dstyaxis);

FX = TriScatteredInterp(dstx, dsty, reshape(srcx,size(dstx)));
FY = TriScatteredInterp(dstx, dsty, reshape(srcy,size(dsty)));

xc = FX(rxp, ryp);
yc = FY(rxp, ryp);


im2 = interp2(srcxaxis, srcyaxis, double(src), xc, yc, '*linear');

imout = reshape(im2, [length(dstyaxis) length(dstxaxis)]);

