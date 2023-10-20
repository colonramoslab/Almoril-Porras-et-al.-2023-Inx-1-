function realpts = realPtsFromCamPts(cc, campts)
%function realpts = realPtsFromCamPts(cc, campts)

x = cc.c2rX(campts(1,:), campts(2,:));
y = cc.c2rY(campts(1,:), campts(2,:));
realpts = [x;y];