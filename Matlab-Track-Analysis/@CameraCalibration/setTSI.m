function setTSI(cc)
%function setTSI(cc)

cc.camx = makecolumn(cc.camx);
cc.camy = makecolumn(cc.camy);
cc.realx = makecolumn(cc.realx);
cc.realy = makecolumn(cc.realy);

[cx, cy, rx, ry] = guessOutsideHull (cc.camx, cc.camy, cc.realx, cc.realy, [min(cc.realx) max(cc.realx)], [min(cc.realy) max(cc.realy)]);
cc.r2cX = TriScatteredInterp (rx, ry, cx);
cc.r2cY = TriScatteredInterp (rx, ry, cy);

[rx, ry, cx, cy] = guessOutsideHull (cc.realx, cc.realy, cc.camx, cc.camy, [min(cc.camx) max(cc.camx)], [min(cc.camy) max(cc.camy)]);
cc.c2rX = TriScatteredInterp (cx, cy, rx);
cc.c2rY = TriScatteredInterp (cx, cy, ry);
%cc.c2rX = TriScatteredInterp (cc.camx, cc.camy, cc.realx);
%cc.c2rY = TriScatteredInterp (cc.camx, cc.camy, cc.realy);

