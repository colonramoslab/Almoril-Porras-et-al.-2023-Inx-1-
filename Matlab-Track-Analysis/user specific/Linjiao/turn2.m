% calculate the error bar of the turning rates by average them over
% experiments (04/26/11)

omega_in=[];
omega_out=[];
omega_early_later=[];
reversal_in=[];
reversal_out=[];
reversal_early_later=[];
turn_delta=[];

% tAxis=150:300:2850;
% fAxis=300:600:5700;

tAxis=60:120:840;
fAxis=120:240:1680;

ne=length(eset.expt);
nb=length(tAxis);
t_turns=cell(1,nb);
t_omega=cell(1,nb);
t_reversal=cell(1,nb);

runlength=[];
runlength_early_later=[];
run_dir=[];
run_long_dir=[];
delta_run=[];
delta_run_early_later=[];


run_up=[];
run_down=[];
run_else=[];





for i=1:length(eset.expt)
     turns=[]; 
     omega=[];
     reversal=[];
     time=eset.expt(i).gatherField('eti');
     for j=1:length(eset.expt(i).track)
         sFrame=eset.expt(i).track(j).startFrame;
         for k=1:length(eset.expt(i).track(j).reorientation)
             reo=eset.expt(i).track(j).reorientation(k);
             reoseq=reo.turnsequence;
             turnf=sFrame+reo.startInd;
             turn_delta=[turn_delta, reo.dTheta];
             for m=1:length(reo.sharpTurn)
                 turns = [turns, (reo.sharpTurn(m).startInd+sFrame)];
                 if reo.sharpTurn(m).typeCode==-1
                     omega = [omega, (reo.sharpTurn(m).startInd+sFrame)];
                 elseif reo.sharpTurn(m).typeCode==1
                     reversal = [reversal, (reo.sharpTurn(m).startInd+sFrame)];
                 end
             end
             if isequal(reoseq, -1)
                 omega_in = [omega_in, reo.thetaIn];
                 omega_out = [omega_out, reo.thetaOut];
                 if turnf<1200
                     omega_early_later=[omega_early_later, 1];
                 elseif turnf>2400 && turnf<3600
                     omega_early_later=[omega_early_later, -1];
                 else
                     omega_early_later=[omega_early_later, 0];
                 end
             elseif isequal(reoseq, [1 -1])
                 reversal_in = [reversal_in, reo.thetaIn];
                 reversal_out = [reversal_out, reo.thetaOut];
                 if turnf<1200
                     reversal_early_later=[reversal_early_later, 1];
                 elseif turnf>2400 && turnf<3600
                     reversal_early_later=[reversal_early_later, -1];
                 else
                     reversal_early_later=[reversal_early_later, 0];
                 end
             end
             runlen=eset.expt(i).track(j).run(k).runTime;
             rundir=eset.expt(i).track(j).run(k).meanTheta;
             runstart=eset.expt(i).track(j).run(k).startInd; 
             runend=eset.expt(i).track(j).run(k).endInd;
             runmid=(runstart+runend)/2;
             runmid=fix(runmid);
             runf=sFrame+runmid;
             if runf<1800 %&& isequal(reoseq, [1 -1])
                 runlength=[runlength, runlen];
             end
             if runf<1200
                     runlength_early_later=[runlength_early_later, 1];
             elseif runf>2400 && runf<3600
                     runlength_early_later=[runlength_early_later, -1];
             else
                     runlength_early_later=[runlength_early_later, 0];
             end

             if runlen>50
                 tracktheta=eset.expt(i).track(j).getDerivedQuantity('theta');
                 runtheta=tracktheta(runmid-10:runmid+10);
                 d_run=runtheta(end)-runtheta(1);
                 delta_run=[delta_run, d_run];
                 if runf<1200
                     delta_run_early_later=[delta_run_early_later, 1];
                 elseif runf>2400 && runf<3600
                     delta_run_early_later=[delta_run_early_later, -1];
                 else
                     delta_run_early_later=[delta_run_early_later, 0];
                 end
                 runtheta=unwrap(runtheta);
                 rundir_long=mean(runtheta);
                 while rundir_long>pi
                     rundir_long=rundir_long-2*pi;
                 end
                 while rundir_long<-pi
                     rundir_long=rundir_long+2*pi;
                 end
                 run_long_dir=[run_long_dir, rundir_long];
             end
             run_dir=[run_dir, rundir];
%              if rundir>-pi/4 && rundir<pi/4 && runf<1800
%                  run_up=[run_up, runlen];
%              elseif rundir<-pi*3/4 || rundir>pi*3/4 && runf<1800
%                  run_down=[run_down, runlen];
%              else
%                  run_else=[run_else, runlen];
%              end
             if rundir>-pi/8 && rundir<pi/8 && runf<1800
                 run_up=[run_up, runlen];
             elseif (rundir<-pi*7/8 || rundir>pi*7/8) && runf<1800
                 run_down=[run_down, runlen];
             elseif ((rundir>pi*3/8 && rundir<pi*5/8) || (rundir>-pi*5/8 && rundir<-pi*3/8)) && runf<1800
                 run_else=[run_else, runlen];
             end
         end
     end

    nt=hist(turns,fAxis);
    nti=hist(time,tAxis);
    turnRate=nt./nti*2; % *2 is to fix the unit from per frame to per sec
    nr=hist(reversal,fAxis);
    revRate=nr./nti*2;
    no=hist(omega,fAxis);
    omegaRate=no./nti*2;
    for n=1:nb
        if nti(n)>600
            t_turns{n}=[t_turns{n}, turnRate(n)];
            t_reversal{n}=[t_reversal{n}, revRate(n)];
            t_omega{n}=[t_omega{n}, omegaRate(n)];
        end
    end
    
end


tRate=cellfun(@mean,t_turns);
e_tRate=cellfun(@std,t_turns)./sqrt(cellfun(@length,t_turns));
figure;
bar(tAxis, tRate);
hold on;
errorbar(tAxis, tRate, e_tRate, '.r');
xlabel ('time(s)');
ylabel('turning rate (per sec)');
title ('Total turning rate vs. time');

rRate=cellfun(@mean,t_reversal);
e_rRate=cellfun(@std,t_reversal)./sqrt(cellfun(@length,t_reversal));
figure;
bar(tAxis, rRate);
hold on;
errorbar(tAxis, rRate, e_rRate, '.r');
xlabel ('time(s)');
ylabel('reversal rate (per sec)');
title ('Reversal rate vs. time');

oRate=cellfun(@mean,t_omega);
e_oRate=cellfun(@std,t_omega)./sqrt(cellfun(@length,t_omega));
figure;
bar(tAxis, oRate);
hold on;
errorbar(tAxis, oRate, e_oRate, '.r');
xlabel ('time(s)');
ylabel('omega turn rate (per sec)');
title ('Omega turn rate vs. time');

figure;
omega_delta=omega_out-omega_in;
for i=1:length(omega_delta)
    if omega_delta(i)>pi
        omega_delta(i)=omega_delta(i)-2*pi;
    elseif omega_delta(i)<-pi
        omega_delta(i)=omega_delta(i)+2*pi;
    end
end
rose(omega_delta,12);
title('one omega turn');

figure;
reversal_delta=reversal_out-reversal_in;
for i=1:length(reversal_delta)
    if reversal_delta(i)>pi
        reversal_delta(i)=reversal_delta(i)-2*pi;
    elseif reversal_delta(i)<-pi
        reversal_delta(i)=reversal_delta(i)+2*pi;
    end
end
rose(reversal_delta,12);
title('reversal omega turn');

omega_early=omega_delta(omega_early_later==1);
omega_later=omega_delta(omega_early_later==-1);
rev_early=reversal_delta(reversal_early_later==1);
rev_later=reversal_delta(reversal_early_later==-1);

omega_delta_1=omega_delta(omega_in<-pi*3/4|omega_in>pi*3/4);
omega_delta_2=omega_delta(omega_in<pi*3/4 & omega_in>pi/4);
omega_delta_3=omega_delta(omega_in<pi/4 & omega_in>-pi/4);
omega_delta_4=omega_delta(omega_in<-pi/4 & omega_in>-pi*3/4);

rev_delta_1=reversal_delta(reversal_in<-pi*3/4|reversal_in>pi*3/4);
rev_delta_2=reversal_delta(reversal_in<pi*3/4 & reversal_in>pi/4);
rev_delta_3=reversal_delta(reversal_in<pi/4 & reversal_in>-pi/4);
rev_delta_4=reversal_delta(reversal_in<-pi/4 & reversal_in>-pi*3/4);

figure;
rose(omega_delta_1,12);
title('hist of delta theta(omega turn, thetaIn<-135 or >135)');
mo1=mean(omega_delta_1);
no1=length(omega_delta_1);
so1=std(omega_delta_1);

figure;
rose(omega_delta_2,12);
title('hist of delta theta(omega turn, thetaIn 45~135)');
mo2=mean(omega_delta_2);
no2=length(omega_delta_2);
so2=std(omega_delta_2);

figure;
rose(omega_delta_3,12);
title('hist of delta theta(omega turn, thetaIn -45~45)');
mo3=mean(omega_delta_3);
no3=length(omega_delta_3);
so3=std(omega_delta_3);

figure;
rose(omega_delta_4,12);
title('hist of delta theta(omega turn, thetaIn -135~-45)');
mo4=mean(omega_delta_4);
no4=length(omega_delta_4);
so4=std(omega_delta_4);

figure;
rose(rev_delta_1,12);
title('hist of delta theta(reversal omega turn, thetaIn<-135 or >135)');
mr1=mean(rev_delta_1);
nr1=length(rev_delta_1);
sr1=std(rev_delta_1);

figure;
rose(rev_delta_2,12);
title('hist of delta theta(reversal omega turn, thetaIn 45~135)');
mr2=mean(rev_delta_2);
nr2=length(rev_delta_2);
sr2=std(rev_delta_2);


figure;
rose(rev_delta_3,12);
title('hist of delta theta(reversal omega turn, thetaIn -45~45)');
mr3=mean(rev_delta_3);
nr3=length(rev_delta_3);
sr3=std(rev_delta_3);


figure;
rose(rev_delta_4,12);
title('hist of delta theta(reversal omega turn, thetaIn -135~-45)');
mr4=mean(rev_delta_4);
nr4=length(rev_delta_4);
sr4=std(rev_delta_4);

figure;
runl=runlength;
runlength=runlength(runlength<100);
runlength_early_later=runlength_early_later(runlength<100);
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
run_up=run_up(run_up<100);
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
run_down=run_down(run_down<100);
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

figure;
run_else=run_else(run_else<100);
[ye,xe]=hist(run_else,tA);
ye = ye/sum(ye);
indx=find(ye>0.01);
xe1=xe(indx)';
ye1=ye(indx)';
xe2=xe1;
ye2=log(ye1);
semilogy(xe1,ye1,'b.','MarkerSize',12);
hold on;
mnxe=min(xe2);
xe2=xe2-mnxe;
a = polyfit(xe2,ye2,1);
fitelse=exp(polyval(a,xe2));
semilogy(xe2+mnxe,fitelse);
xlabel ('Duration(s)');
ylabel ('Count(%)');
title (['hist of runlength orthogonal, tau =' num2str(-1/a(1))]);

% figure;
% rose(delta_run,12)
% title('hist of delta theta within one run');
% 
% for i=1:length(delta_run)
%     if delta_run(i)>pi
%         delta_run(i)=delta_run(i)-2*pi;
%     elseif delta_run(i)<-pi
%         delta_run(i)=delta_run(i)+2*pi;
%     end
% end
% 
% delta_run1=delta_run(run_long_dir<-pi*3/4|run_long_dir>pi*3/4);
% delta_run2=delta_run(run_long_dir<pi*3/4 & run_long_dir>pi/4);
% delta_run3=delta_run(run_long_dir<pi/4 & run_long_dir>-pi/4);
% delta_run4=delta_run(run_long_dir<-pi/4 & run_long_dir>-pi*3/4);
% 
% figure;
% rose(delta_run1,12);
% title('hist of delta theta within one run(<-135 or >135)');
% m1=mean(delta_run1);
% n1=length(delta_run1);
% s1=std(delta_run1);
% 
% figure;
% rose(delta_run2,12);
% title('hist of delta theta within one run(45~135)');
% m2=mean(delta_run2);
% n2=length(delta_run2);
% s2=std(delta_run2);
% 
% figure;
% rose(delta_run3,12);
% title('hist of delta theta within one run(-45~45)');
% m3=mean(delta_run3);
% n3=length(delta_run3);
% s3=std(delta_run3);
% 
% figure;
% rose(delta_run4,12);
% title('hist of delta theta within one run(-135~-45)');
% m4=mean(delta_run4);
% n4=length(delta_run4);
% s4=std(delta_run4);


             
         