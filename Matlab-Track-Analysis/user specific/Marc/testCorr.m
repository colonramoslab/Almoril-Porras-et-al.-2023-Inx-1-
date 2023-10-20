
vfield = 'nvel';
hfield = 'smhdir';
tfield = 'stmdir';
figure(1);
[xc_wr, np_wr, tx_wr] = e1.crosscorrelate(vfield, hfield, 'withinRuns', true);
[xc_ir, np_ir, tx_ir] = e1.crosscorrelate(vfield, hfield, 'inRuns', true);
[xc_h, np_h, tx_h] = e1.crosscorrelate(vfield, hfield);

semilogy(tx_wr, xc_wr./np_wr, 'b.-', tx_ir, xc_ir./np_ir, 'g.-', tx_h, xc_h./np_h, 'r.-'); axis([-5 5 .9 1]);
legend('within runs', 'in runs', 'normal');

figure(2);
[xc_h, np_h, tx_h] = e1.crosscorrelate(vfield, hfield, 'withinRuns', true);
[xc_t, np_t, tx_t] = e1.crosscorrelate(vfield, tfield, 'withinRuns', true);
[xc_ht, np_ht, tx_ht] =  e1.crosscorrelate(tfield, hfield, 'withinRuns', true);
semilogy(tx_h, xc_h./np_h, 'g.-', tx_t, xc_t./np_t, 'r.-', tx_ht, xc_ht./np_ht, 'b.-'); axis([-15 15 .9 1]);
legend('vel-head', 'vel-tail', 'tail-head');
title ('cross correlation within runs');

figure(3);
[ac_v, np_v, tx_v] = e1.autocorrelate(vfield, 'withinRuns', true);
[ac_h, np_h, tx_h] = e1.autocorrelate(hfield, 'withinRuns', true);
[ac_t, np_t, tx_t] = e1.autocorrelate(tfield, 'withinRuns', true);
semilogy(tx_v, ac_v./np_v, 'b.', tx_h, ac_h./np_h, 'g.', tx_t, ac_t./np_t, 'r.'); axis([0 20 0.8 1]);
title ('auto correlation within runs'); xlabel('tau(s)');
legend ('vel', 'head dir', 'tail dir');
return




nsteps = 1E3;
ncircles = 1;
t = (0:1:nsteps)*ncircles*2*pi/nsteps;
x = cos(t);
y = sin(t);
loc = [x;y];
size(loc)
sloc = lowpass1d(loc, 10);
vel = deriv(loc, 1);
svel = deriv(sloc, 1);
l = sqrt(sum(vel.^2));
vnorm = vel./[l;l];
l = sqrt(sum(svel.^2));
svnorm = svel./[l;l];

%[cunb,c, np] = xcorrVec (a, b)
figure(1);
[cunb,c,np] = xcorrVec(vnorm,svnorm);
tx = (1:length(cunb))-(length(cunb)+1)/2;
semilogy (tx, cunb,'b.-'); xlim([-100, 100]); ylim([0.9 1]);

headang = t+pi/2 + pi/18*randn(size(t));
headdir = [cos(headang);sin(headang)];
bl = 0.2;
hloc = loc + bl*headdir;
shead = lowpass1d(hloc, 10);

l = sqrt(sum(vel.^2));
vnorm = vel./[l;l];

figure(2);
inds = 1:25:nsteps;
plot (loc(1,:), loc(2,:), 'b-');%, hloc(1,:), hloc(2,:), 'r-')
hold on
quiver(loc(1,inds), loc(2,inds), shead(1,inds)-sloc(1,inds), shead(2,inds)-sloc(2,inds),0,'r');
quiver(loc(1,inds), loc(2,inds), vnorm(1,inds), vnorm(2,inds),0,'y');

hold off
axis equal

figure(3);
mh = shead-sloc;
l = sqrt(sum(mh.^2));
mhdir = mh./[l;l];
[cunb,c,np] = xcorrVec(vnorm,mhdir);
tx = (1:length(cunb))-(length(cunb)+1)/2;
semilogy (tx, cunb,'b.-'); xlim([-10, 10]); ylim([0.99 1]); 
