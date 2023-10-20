function magfactor = pixelsPerRealUnit (cc)
%function magfactor = pixelsPerRealUnit (cc)
%
% how many pixels = 1 real unit (e.g. cm)
magfactor = 1/cc.realUnitsPerPixel;
