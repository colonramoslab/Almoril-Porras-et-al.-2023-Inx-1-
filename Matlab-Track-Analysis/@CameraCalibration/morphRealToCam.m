function camim = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin)
%function camin = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin)
%
%optional args 'camxaxis', 'camyaxis'
%whos cc
camxaxis = min(cc.camx):max(cc.camx);
camyaxis = min(cc.camy):max(cc.camy);
%camxaxis = 1:size(camim,2);
%camyaxis = 1:size(camim,1);

varargin = assignApplicable(varargin);

    
[cxp,cyp] = meshgrid(camxaxis,camyaxis);
cxp = cxp(:);
cyp = cyp(:);

xr = cc.c2rX(cxp, cyp);
yr = cc.c2rY(cxp, cyp);

im2 = interp2(realxaxis, realyaxis, double(realim), xr, yr, '*linear');

camim = reshape(im2, [length(camyaxis) length(camxaxis)]);
