%calibrating etac ethyl acetate triangle to produce linear ramp from gas sensor

%data = importdata('\\LABNAS2\LarvalCO2\gas calibrations\etac triangle wave\inlet trace of etac with pure substance and ten minute triangle wave good detector.txt');
data = importdata('\\LABNAS2\LarvalCO2\gas calibrations\etac triangle wave\inlet trace of etac with pure substance and ten minute triangle with pure substance using triangle_target_rev_1.txt');

pidt = data.data(:,2);
pidppm = data.data(:,1);
gast = data.data(:,4);
gassp = data.data(:,3);

%t0 = max(min(pidt), min(gast));
%t1 = min(max(pidt), max(gast));
%t0 = 680E3;
%t1 = t0 + 40*60E3;
t0 = 330E3;
t1 = t0 + 30*60E3;

dt = 35;%min(diff(unique(pidt(isfinite(pidt)))));
tx = t0:dt:t1;

inds = isfinite(pidt) & isfinite(pidppm);
pidt = pidt(inds);
pidppm = pidppm(inds);
[pidt,I] = unique(pidt);
pidppm = pidppm(I);

inds = isfinite(gast) & isfinite(gassp);
gast = gast(inds);
gassp = gassp(inds);
[gast,I] = unique(gast);
gassp = gassp(I);

ppm = lowpass1D(interp1(pidt, pidppm, tx),10);
sp = interp1(gast, gassp, tx);

figure(1); clf;
plot (tx, sp/max(sp), tx, ppm/max(ppm));

figure(2); clf
[c,lags] = xcorr((sp - mean(sp))/max(sp), (ppm-mean(ppm))/max(ppm), ceil(90000/dt), 'unbiased');
plot (lags*dt/1000, c);

[~,I] = max(c);
delay_t = dt*lags(I);

figure(3); clf;
ppmshift = interp1(tx+delay_t, ppm, tx);
plot (tx, sp/max(sp), tx, ppmshift/max(ppm));

return;

figure(4); clf;
%plot (sp([true diff(sp)>0]), ppmshift([true diff(sp)>0]), 'r-', sp([false diff(sp)<0]), ppmshift([false diff(sp)<0]), 'b-');

ppmup = ppmshift(diff(sp) > 0);
spup = sp(diff(sp) > 0);
[spup,I] = sort(spup);
ppmup = ppmup(I);

ppmdn = ppmshift(diff(sp) < 0);
spdn = sp(diff(sp) < 0);
[spdn,I] = sort(spdn);
ppmdn = ppmdn(I);
spx = 0:0.1:50;
[spu, pu] = meanyvsx(spup, ppmup, spx);
[spd, pd] = meanyvsx(spdn, ppmdn, spx);

inds = isfinite(spu) & isfinite(pu);
pu = pu(inds);
spu = spu(inds);

inds = isfinite(spd) & isfinite(pd);
pd = pd(inds);
spd = spd(inds);



plot (spup, ppmup, 'k.', spu, pu, 'r-', spdn, ppmdn,'g.',spd,pd,'b-', 'LineWidth', 4, 'MarkerSize', 1);

maxp = min(max(pu), max(pd));
minp = max(min(pu), min(pd));
sptarget_up = interp1(pu,spu, linspace(minp, maxp, 300));
sptarget_dn = interp1(pd,spd, linspace(maxp, minp, 300));

sptarget = [sptarget_up sptarget_dn];
sp_tx = (1:length(sptarget))- 1;


figure(5); clf
plot (sp_tx,sptarget);



%csvwrite('\\LABNAS2\LarvalCO2\gas calibrations\etac triangle wave\triangle_target_rev_1.txt', [sp_tx;sptarget]');