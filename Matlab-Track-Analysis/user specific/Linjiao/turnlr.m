runl1=delta_run1(delta_run1>0); 
runr1=-delta_run1(delta_run1<0);
runl2=delta_run2(delta_run2>0); 
runr2=-delta_run2(delta_run2<0);
runl3=delta_run3(delta_run3>0); 
runr3=-delta_run3(delta_run3<0);
runl4=delta_run4(delta_run4>0); 
runr4=-delta_run4(delta_run4<0);

prun1=length(runl1)/length(delta_run1);
prun2=length(runl2)/length(delta_run2);
prun3=length(runl3)/length(delta_run3);
prun4=length(runl4)/length(delta_run4);
erun1=sqrt(prun1*(1-prun1)/length(delta_run1));
erun2=sqrt(prun2*(1-prun2)/length(delta_run2));
erun3=sqrt(prun3*(1-prun3)/length(delta_run3));
erun4=sqrt(prun4*(1-prun4)/length(delta_run4));

prun=[prun1 prun2 prun3 prun4];
erun=[erun1 erun2 erun3 erun4];
figure;
bar(prun);
hold on
errorbar(prun, erun, 'xr');
title('probility of left run steering');

mrunl1=mean(runl1);
mrunr1=mean(runr1);
mrunl2=mean(runl2);
mrunr2=mean(runr2);
mrunl3=mean(runl3);
mrunr3=mean(runr3);
mrunl4=mean(runl4);
mrunr4=mean(runr4);

erunl1=std(runl1)/sqrt(length(runl1));
erunr1=std(runr1)/sqrt(length(runr1));
erunl2=std(runl2)/sqrt(length(runl2));
erunr2=std(runr2)/sqrt(length(runr2));
erunl3=std(runl3)/sqrt(length(runl3));
erunr3=std(runr3)/sqrt(length(runr3));
erunl4=std(runl4)/sqrt(length(runl4));
erunr4=std(runr4)/sqrt(length(runr4));

mrun=[mrunl1 mrunr1;mrunl2 mrunr2;mrunl3 mrunr3;mrunl4 mrunr4];
erun=[erunl1 erunr1;erunl2 erunr2;erunl3 erunr3;erunl4 erunr4];
mrun=rad2deg(mrun);
erun=rad2deg(erun);
figure;
groupederrorbar(mrun, erun);
legend('left','right');
title('mean of left right run steering');



ol1=omega_delta_1(omega_delta_1>0); 
or1=-omega_delta_1(omega_delta_1<0);
ol2=omega_delta_2(omega_delta_2>0); 
or2=-omega_delta_2(omega_delta_2<0);
ol3=omega_delta_3(omega_delta_3>0); 
or3=-omega_delta_3(omega_delta_3<0);
ol4=omega_delta_4(omega_delta_4>0); 
or4=-omega_delta_4(omega_delta_4<0);

rl1=rev_delta_1(rev_delta_1>0); 
rr1=-rev_delta_1(rev_delta_1<0);
rl2=rev_delta_2(rev_delta_2>0); 
rr2=-rev_delta_2(rev_delta_2<0);
rl3=rev_delta_3(rev_delta_3>0); 
rr3=-rev_delta_3(rev_delta_3<0);
rl4=rev_delta_4(rev_delta_4>0); 
rr4=-rev_delta_4(rev_delta_4<0);

po1=length(ol1)/length(omega_delta_1);
po2=length(ol2)/length(omega_delta_2);
po3=length(ol3)/length(omega_delta_3);
po4=length(ol4)/length(omega_delta_4);
eo1=sqrt(po1*(1-po1)/length(omega_delta_1));
eo2=sqrt(po2*(1-po2)/length(omega_delta_2));
eo3=sqrt(po3*(1-po3)/length(omega_delta_3));
eo4=sqrt(po4*(1-po4)/length(omega_delta_4));

no1=length(omega_delta_1);
ko1=max(length(ol1),length(or1));
no2=length(omega_delta_2);
ko2=max(length(ol2),length(or2));
no3=length(omega_delta_3);
ko3=max(length(ol3),length(or3));
no4=length(omega_delta_4);
ko4=max(length(ol4),length(or4));

nr1=length(rev_delta_1);
kr1=max(length(rl1),length(rr1));
nr2=length(rev_delta_2);
kr2=max(length(rl2),length(rr2));
nr3=length(rev_delta_3);
kr3=max(length(rl3),length(rr3));
nr4=length(rev_delta_4);
kr4=max(length(rl4),length(rr4));

po=[po1 po2 po3 po4];
eo=[eo1 eo2 eo3 eo4];
figure;
bar(po);
hold on
errorbar(po, eo, 'xr');
title('probility of left turns (omega)');

pr1=length(rl1)/length(rev_delta_1);
pr2=length(rl2)/length(rev_delta_2);
pr3=length(rl3)/length(rev_delta_3);
pr4=length(rl4)/length(rev_delta_4);
er1=sqrt(pr1*(1-pr1)/length(rev_delta_1));
er2=sqrt(pr2*(1-pr2)/length(rev_delta_2));
er3=sqrt(pr3*(1-pr3)/length(rev_delta_3));
er4=sqrt(pr4*(1-pr4)/length(rev_delta_4));

pr=[pr1 pr2 pr3 pr4];
er=[er1 er2 er3 er4];
figure;
bar(pr);
hold on
errorbar(pr, er, 'xr');
title('probility of left turns (reversal omega)');

mol1=mean(ol1);
sol1=std(ol1);
nol1=length(ol1);

mor1=mean(or1);
sor1=std(or1);
nor1=length(or1);

mol2=mean(ol2);
sol2=std(ol2);
nol2=length(ol2);

mor2=mean(or2);
sor2=std(or2);
nor2=length(or2);

mol3=mean(ol3);
sol3=std(ol3);
nol3=length(ol3);

mor3=mean(or3);
sor3=std(or3);
nor3=length(or3);

mol4=mean(ol4);
sol4=std(ol4);
nol4=length(ol4);

mor4=mean(or4);
sor4=std(or4);
nor4=length(or4);

mrl1=mean(rl1);
srl1=std(rl1);
nrl1=length(rl1);

mrr1=mean(rr1);
srr1=std(rr1);
nrr1=length(rr1);

mrl2=mean(rl2);
srl2=std(rl2);
nrl2=length(rl2);

mrr2=mean(rr2);
srr2=std(rr2);
nrr2=length(rr2);

mrl3=mean(rl3);
srl3=std(rl3);
nrl3=length(rl3);

mrr3=mean(rr3);
srr3=std(rr3);
nrr3=length(rr3);

mrl4=mean(rl4);
srl4=std(rl4);
nrl4=length(rl4);

mrr4=mean(rr4);
srr4=std(rr4);
nrr4=length(rr4);


eol1=std(ol1)/sqrt(length(ol1));
eor1=std(or1)/sqrt(length(or1));
eol2=std(ol2)/sqrt(length(ol2));
eor2=std(or2)/sqrt(length(or2));
eol3=std(ol3)/sqrt(length(ol3));
eor3=std(or3)/sqrt(length(or3));
eol4=std(ol4)/sqrt(length(ol4));
eor4=std(or4)/sqrt(length(or4));

erl1=std(rl1)/sqrt(length(rl1));
err1=std(rr1)/sqrt(length(rr1));
erl2=std(rl2)/sqrt(length(rl2));
err2=std(rr2)/sqrt(length(rr2));
erl3=std(rl3)/sqrt(length(rl3));
err3=std(rr3)/sqrt(length(rr3));
erl4=std(rl4)/sqrt(length(rl4));
err4=std(rr4)/sqrt(length(rr4));

mo=[mol1 mor1;mol2 mor2;mol3 mor3;mol4 mor4];
eo=[eol1 eor1;eol2 eor2;eol3 eor3;eol4 eor4];
mo=rad2deg(mo);
eo=rad2deg(eo);
figure;
groupederrorbar(mo, eo);
legend('left','right');
title('mean of left right turns(omega)');

mr=[mrl1 mrr1;mrl2 mrr2;mrl3 mrr3;mrl4 mrr4];
er=[erl1 err1;erl2 err2;erl3 err3;erl4 err4];
mr=rad2deg(mr);
er=rad2deg(er);
figure;
groupederrorbar(mr, er);
legend('left','right');
title('mean of left right turns(reversal omega)');

