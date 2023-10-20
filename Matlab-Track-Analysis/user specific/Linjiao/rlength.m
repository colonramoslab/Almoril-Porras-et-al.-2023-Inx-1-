

runlength=[];
run_up=[];
run_down=[];
run_else=[];





for i=1:length(eset.expt)
     for j=1:length(eset.expt(i).track)
         sFrame=eset.expt(i).track(j).startFrame;
         if eset.expt(i).track(j).npts>3600 && eset.expt(i).track(j).startFrame<120 % add 06/28/11
         for k=1:length(eset.expt(i).track(j).reorientation)
             runlen=eset.expt(i).track(j).run(k).runTime;
             rundir=eset.expt(i).track(j).run(k).meanTheta;
             runlength=[runlength, runlen];
             runstart=eset.expt(i).track(j).run(k).startInd; 
             runend=eset.expt(i).track(j).run(k).endInd;
             runmid=(runstart+runend)/2;
             runmid=fix(runmid);
             runf=sFrame+runmid;
             run_dir=[run_dir, rundir];
             if rundir>-pi/4 && rundir<pi/4 && runf<1800
                 run_up=[run_up, runlen];
             elseif rundir<-pi*3/4 || rundir>pi*3/4 && runf<1800
                 run_down=[run_down, runlen];
             elseif runf<1800
                 run_else=[run_else, runlen];
             end
         end
         end
     end
    
end




figure;
runlength=runlength(runlength<600);
tA=7.5:5:597.5;
[y,x]=hist(runlength,tA);
y = y/sum(y);
indx=find(y>0.01);
x1=x(indx)';
y1=y(indx)';
x2=x1;
y2=log(y1);
semilogy(x1,y1,'b.','MarkerSize',12);
hold on;
mnx=min(x2);
x2=x2-mnx;
a = polyfit(x2,y2,1);
fit=exp(polyval(a,x2));
semilogy(x2+mnx,fit);
xlabel ('Duration(s)');
ylabel ('Count(%)');
title (['hist of runlength, tau =' num2str(-1/a(1))]);

figure;
run_up=run_up(run_up<600);
[yu,xu]=hist(run_up,tA);
yu = yu/sum(yu);
indx=find(yu>0.01);
xu1=xu(indx)';
yu1=yu(indx)';
xu2=xu1;
yu2=log(yu1);
semilogy(xu1,yu1,'b.','MarkerSize',12);
hold on;
mnxu=min(xu2);
xu2=xu2-mnxu;
a = polyfit(xu2,yu2,1);
fitu=exp(polyval(a,xu2));
semilogy(xu2+mnxu,fitu);
xlabel ('Duration(s)');
ylabel ('Count(%)');
title (['hist of runlength up gradient, tau =' num2str(-1/a(1))]);

figure;
run_down=run_down(run_down<600);
[yd,xd]=hist(run_down,tA);
yd = yd/sum(yd);
indx=find(yd>0.01);
xd1=xd(indx)';
yd1=yd(indx)';
xd2=xd1;
yd2=log(yd1);
semilogy(xd1,yd1,'b.','MarkerSize',12);
hold on;
mnxd=min(xd2);
xd2=xd2-mnxd;
a = polyfit(xd2,yd2,1);
fitdown=exp(polyval(a,xd2));
semilogy(xd2+mnxd,fitdown);
xlabel ('Duration(s)');
ylabel ('Count(%)');
title (['hist of runlength down gradient, tau =' num2str(-1/a(1))]);

