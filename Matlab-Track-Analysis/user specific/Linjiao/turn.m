turns=[]; 
omega=[];
reversal=[];
omega_in=[];
omega_out=[];
reversal_in=[];
reversal_out=[];
time=eset.gatherField('eti');

for i=1:length(eset.expt)
     for j=1:length(eset.expt(i).track)
         for k=1:length(eset.expt(i).track(j).reorientation)
             reo=eset.expt(i).track(j).reorientation(k);
             reoseq=reo.turnsequence;
             for m=1:length(reo.sharpTurn)
                 turns = [turns, reo.sharpTurn(m).startInd];
                 if reo.sharpTurn(m).typeCode==-1
                     omega = [omega, reo.sharpTurn(m).startInd];
                 elseif reo.sharpTurn(m).typeCode==1
                     reversal = [reversal, reo.sharpTurn(m).startInd];
                 end
             end
             if isequal(reoseq, -1)
                 omega_in = [omega_in, reo.thetaIn];
                 omega_out = [omega_out, reo.thetaOut];
             elseif isequal(reoseq, [1 -1])
                 reversal_in = [reversal_in, reo.thetaIn];
                 reversal_out = [reversal_out, reo.thetaOut];
             end
         end
     end
end

tAxis=150:300:2850;
fAxis=300:600:5700;
nt=hist(turns,fAxis);
nti=hist(time,tAxis);
turnRate=nt./nti*2; % *2 is to fix the unit from per frame to per sec
figure;
bar(tAxis, turnRate);
xlabel ('time(s)');
ylabel('turning rate (per sec)');
title ('Total turning rate vs. time');

nr=hist(reversal,fAxis);
revRate=nr./nti*2;
figure;
bar(tAxis, revRate);
xlabel ('time(s)');
ylabel('reversal rate (per sec)');
title ('Reversal rate vs. time');

no=hist(omega,fAxis);
omegaRate=no./nti*2;
figure;
bar(tAxis, omegaRate);
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
omega_delta=rad2deg(omega_delta);
omega_in=rad2deg(omega_in);
plot(omega_delta,omega_in,'.','MarkerSize',10);
xlabel ('delta theta)');
ylabel('theta in');
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
reversal_delta=rad2deg(reversal_delta);
reversal_in=rad2deg(reversal_in);
plot(reversal_delta,reversal_in,'.','MarkerSize',10);
xlabel ('delta theta');
ylabel('theta in)');
title('reversal omega turn');

omega_delta_1=omega_delta(find(omega_in<-135|omega_in>135));
omega_delta_2=omega_delta(find(omega_in<135 & omega_in>45));
omega_delta_3=omega_delta(find(omega_in<45 & omega_in>-45));
omega_delta_4=omega_delta(find(omega_in<-45 & omega_in>-135));

rev_delta_1=reversal_delta(find(reversal_in<-135|reversal_in>135));
rev_delta_2=reversal_delta(find(reversal_in<135 & reversal_in>45));
rev_delta_3=reversal_delta(find(reversal_in<45 & reversal_in>-45));
rev_delta_4=reversal_delta(find(reversal_in<-45 & reversal_in>-135));

xAxis = -165:30:165;
o1=hist(omega_delta_1,xAxis);
o1=o1./sum(o1);
figure;
bar(xAxis,o1);
title('hist of delta theta(omega turn, thetaIn<-135 or >135)');

o2=hist(omega_delta_2,xAxis);
o2=o2./sum(o2);
figure;
bar(xAxis,o2);
title('hist of delta theta(omega turn, thetaIn 45~135)');

o3=hist(omega_delta_3,xAxis);
o3=o3./sum(o3);
figure;
bar(xAxis,o3);
title('hist of delta theta(omega turn, thetaIn -45~45)');

o4=hist(omega_delta_4,xAxis);
o4=o4./sum(o4);
figure;
bar(xAxis,o4);
title('hist of delta theta(omega turn, thetaIn -135~-45)');

r1=hist(rev_delta_1,xAxis);
r1=r1./sum(r1);
figure;
bar(xAxis,r1);
title('hist of delta theta(reversal omega turn, thetaIn<-135 or >135)');

r2=hist(rev_delta_2,xAxis);
r2=r2./sum(r2);
figure;
bar(xAxis,r2);
title('hist of delta theta(reversal omega turn, thetaIn 45~135)');

r3=hist(rev_delta_3,xAxis);
r3=r3./sum(r3);
figure;
bar(xAxis,r3);
title('hist of delta theta(reversal omega turn, thetaIn -45~45)');

r4=hist(rev_delta_4,xAxis);
r4=r4./sum(r4);
figure;
bar(xAxis,r4);
title('hist of delta theta(reversal omega turn, thetaIn -135~-45)');

             
             
         