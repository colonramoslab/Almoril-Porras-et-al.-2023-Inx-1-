function atp = populateMaggotFields (atp)

nptsInSpine = 11;

mcdf = atp.mcdf;
atp.head = mcdf.Head' * AndyTrackPoint.mmPerPixel();
atp.tail = mcdf.Tail' * AndyTrackPoint.mmPerPixel();
atp.et = mcdf.TimeElapsed;
x = [mcdf.BoundaryA(1:2:end) mcdf.BoundaryB(end-1:-2:1)] * AndyTrackPoint.mmPerPixel;
y = [mcdf.BoundaryA(2:2:end) mcdf.BoundaryB(end:-2:2)] * AndyTrackPoint.mmPerPixel;
atp.contour = [x;y];

[geom,iner] = polygeom(x,y);
atp.area = geom(1);
atp.loc = [geom(2);geom(3)];
atp.cov = [iner(4); iner(6); iner(5)]; 

x = mcdf.SegmentedCenterline(1:2:end) * AndyTrackPoint.mmPerPixel;
y = mcdf.SegmentedCenterline(2:2:end) * AndyTrackPoint.mmPerPixel;
d = sqrt((x - x(1)).^2 + (y-y(1)).^2);
[d,I] = unique(d);
x = x(I);
y = y(I);
dx = linspace(max(d), 0, nptsInSpine);
x = interp1(d,x,dx,'linear');
y = interp1(d,y,dx,'linear');
atp.spine = [x;y];
atp.mid = atp.spine(:,ceil(nptsInSpine/2));
atp.htValid = true;
