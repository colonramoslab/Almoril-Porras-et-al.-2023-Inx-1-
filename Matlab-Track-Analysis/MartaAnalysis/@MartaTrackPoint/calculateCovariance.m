function cov = calculateCovariance(mtp)




%%[cc,d] = listOfPointsToBWCC(mtp.contour);
%%rp = regionprops(cc, 'MajorAxisLength','MinorAxisLength','Orientation')
[~,iner] = polygeom(mtp.contour(1,:), mtp.contour(2,:));
%{
covd = [rp.MajorAxisLength*d,0;0,rp.MinorAxisLength*d];
rmat = [cosd(rp.Orientation),sind(rp.Orientation);-sind(rp.Orientation),cosd(rp.Orientation)];
rmat2 = [cosd(rp.Orientation),-sind(rp.Orientation);sind(rp.Orientation),cosd(rp.Orientation)];
cov = rmat*(covd.^2)*rmat2;
cov = [cov(1);cov(2);cov(4)];
%}
cov = [iner(4); iner(6); iner(5)]; 