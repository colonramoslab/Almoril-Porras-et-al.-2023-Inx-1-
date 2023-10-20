gq = eset.expt(2).globalQuantity;

et1 = min(gq(1).xData.et):0.2:max(gq(1).xData.et);
etextra = 0.2:0.2:500;
sp1 = interp1(gq(1).xData.et, gq(1).yData, et1);
spextra = interp1(gq(1).xData.et, gq(1).yData, et1(end) + etextra - 599);
et1 = [et1 etextra+et1(end)];
sp1 = [sp1 spextra];

et1 = et1-min(et1);


plot (et1, sp1); 

vet1 = min(gq(3).xData.et):0.2:(max(gq(3).xData.et) - 30);
vetextra = 0.2:0.2:600;
voc1 = interp1(gq(3).xData.et, gq(3).yData, vet1);
vocextra = interp1(gq(3).xData.et, gq(3).yData, vet1(end) + vetextra - 599);
vet1 = [vet1 vetextra+vet1(end)];
voc1 = [voc1 vocextra];

figure;
plot (vet1, voc1);

gqgas = gq(3);
gqgas.xData.et = vet1;
gqgas.yData = voc1;
for j = 3:6
    gqn = eset.expt(j).globalQuantity;
    etn = min(gqn(1).xData.et):0.2:max(gqn(1).xData.et);
    spn = interp1(gqn(1).xData.et, gqn(1).yData, etn);
    
    etn = etn-min(etn);
    
    inds = 1:min(length(et1), length(etn));
   
    y1 = sp1(inds);
    y2 = spn(inds);
    inds = 50:(inds(end) - 50);
    myfun = @(x) sum((y1(inds) - interp1(y2, inds + x)).^2);
    x = fminsearch(myfun, 0);
    
%    [c,lags] = xcorr((sp1(inds) - mean(sp1))/max(sp1), (spn(inds) - mean(spn))/max(spn), 300, 'unbiased');
 %   [~,I] = max(c);
  %  dt = lags(I)*0.2;
    dt = x/5
    plot (et1(inds)+dt, sp1(inds), etn(inds) , spn(inds)); 
    gqnew = gqgas;
    gqnew.xData.et = gqgas.xData.et + dt;
    eset.expt(j).addGlobalQuantity(gqnew);
end
    