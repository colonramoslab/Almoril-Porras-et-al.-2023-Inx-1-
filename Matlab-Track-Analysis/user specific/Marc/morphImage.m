function imout = morphImage(im, imxaxis, imyaxis, realx, realy, camerax, cameray, realxaxis, realyaxis)
%function imout = morphImage(im, imxaxis, imyaxis, realx, realy, camerax, cameray, realxaxis, realyaxis)
%
%image is a picture taken in camera coordinates
%camerax,cameray <---> realx,realy are corresponding camera points and real
%points
%realxaxis,realyaxis are axes in real coordinates that the image is to be
%morphed to;  note that this may be any subset of the locations spanned by
%realx,realy
if (size(realx,1) == 1)
    realx = realx';
end
if (size(realy,1) == 1)
    realy = realy';
end
if (size(camerax,1) == 1)
    camerax = camerax';
end
if (size(cameray,1) == 1)
    cameray = cameray';
end


if (isempty(imxaxis))
    imxaxis = 1:size(im,2);
end
if (isempty(imyaxis))
    imyaxis = 1:size(im,1);
end

if (isempty(realxaxis))
    realxaxis = linspace(min(realx), max(realx), length(imxaxis));
end
if (isempty(realyaxis))
    realyaxis = linspace(min(realy), max(realy), length(imyaxis));
end


[rxp,ryp] = meshgrid(realxaxis,realyaxis);


imxaxis(1)
imxaxis(2)
imxaxis(end)
imyaxis(1)
imyaxis(2)
imyaxis(end)

FX = TriScatteredInterp(realx, realy, reshape(camerax,size(realx)));
FY = TriScatteredInterp(realx, realy, reshape(cameray,size(realy)));

xc = FX(rxp, ryp);
yc = FY(rxp, ryp);

xc(1)
xc(2)
xc(end)
yc(1)
yc(2)
yc(end)


im2 = interp2(imxaxis, imyaxis, double(im), xc, yc, '*linear');

imout = reshape(im2, [length(realyaxis) length(realxaxis)]);