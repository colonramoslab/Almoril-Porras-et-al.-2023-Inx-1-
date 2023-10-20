coin_flips=rand(100,10000);
coin_flips=coin_flips>0.5;
nHeads=sum(coin_flips,1);
hist(nHeads,0:1:100);
mean(nHeads)
std(nHeads)