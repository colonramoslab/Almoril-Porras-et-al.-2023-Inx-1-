turn_dtheta=[];
turn_raising=[];

omega_dtheta=[];
omega_raising=[];

reversal_dtheta=[];
reversal_raising=[];


for i=1:length(eset.expt)
    time=eset.expt(i).globalQuantity(2).xData;
    dtemp=eset.expt(i).globalQuantity(2).yData;
     for j=1:length(eset.expt(i).track)
         for k=1:length(eset.expt(i).track(j).reorientation)
             reo=eset.expt(i).track(j).reorientation(k);
             reoseq=reo.turnsequence;
             reostart=reo.startInd;
             reoend=reo.endInd;
             tmin=eset.expt(i).track(j).pt(reostart).et;
             tmax=eset.expt(i).track(j).pt(reoend).et;
             fmin=fix(14403/3600.5*tmin);
             fmax=fix(14403/3600.5*tmax);
             turn_dtemp=dtemp(fmin:fmax);
             t_raising=turn_dtemp>0;
             if mean(t_raising)==1 
                t_ind=1; %raising phase
             elseif mean(t_raising)==0
                 t_ind =-1; %falling phase
             else
                 t_ind=0;
             end 
             turn_dtheta=[turn_dtheta, reo.dTheta];
             turn_raising=[turn_raising, t_ind];
             if isequal(reoseq, -1)
                omega_dtheta=[omega_dtheta, reo.dTheta];
                omega_raising=[omega_raising, t_ind];
             elseif isequal(reoseq, [1 -1])
                reversal_dtheta=[reversal_dtheta, reo.dTheta];
                reversal_raising=[reversal_raising, t_ind];
             end
         end
     end
end

figure;
for i=1:length(turn_dtheta)
    if turn_dtheta(i)>pi
        turn_dtheta(i)=turn_dtheta(i)-2*pi;
    elseif turn_dtheta(i)<-pi
        turn_dtheta(i)=turn_dtheta(i)+2*pi;
    end
end
rose(turn_dtheta,12);
title('all turns');

turn_r_dtheta=turn_dtheta(turn_raising==1);
turn_f_dtheta=turn_dtheta(turn_raising==-1);
figure;
rose(turn_r_dtheta,12);
title('all turns, raising temperature');
figure;
rose(turn_f_dtheta,12);
title('all turns, falling temperature');



figure;
for i=1:length(omega_dtheta)
    if omega_dtheta(i)>pi
        omega_dtheta(i)=omega_dtheta(i)-2*pi;
    elseif omega_dtheta(i)<-pi
        omega_dtheta(i)=omega_dtheta(i)+2*pi;
    end
end
rose(omega_dtheta,12);
title('one omega turn');

omega_r_dtheta=omega_dtheta(omega_raising==1);
omega_f_dtheta=omega_dtheta(omega_raising==-1);
figure;
rose(omega_r_dtheta,12);
title('one omega turn, raising temperature');
figure;
rose(omega_f_dtheta,12);
title('one omega turn, falling temperature');

figure;
for i=1:length(reversal_dtheta)
    if reversal_dtheta(i)>pi
        reversal_dtheta(i)=reversal_dtheta(i)-2*pi;
    elseif reversal_dtheta(i)<-pi
        reversal_dtheta(i)=reversal_dtheta(i)+2*pi;
    end
end
rose(reversal_dtheta,12);
title('reversal omega turn');

reversal_r_dtheta=reversal_dtheta(reversal_raising==1);
reversal_f_dtheta=reversal_dtheta(reversal_raising==-1);
figure;
rose(reversal_r_dtheta,12);
title('reversal omega turn, raising temperature');
figure;
rose(reversal_f_dtheta,12);
title('reversal omega turn, falling temperature');

             
         