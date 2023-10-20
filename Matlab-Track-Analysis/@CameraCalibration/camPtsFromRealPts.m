function campts = camPtsFromRealPts(cc, realpts)
%function campts = camPtsFromRealPts(cc, realpts)
realpts = double(realpts);
x = cc.r2cX(realpts(1,:), realpts(2,:));
y = cc.r2cY(realpts(1,:), realpts(2,:));
campts = [x;y];