% calculate averaged run speed in 4 quadrants

run1=[];
run2=[];
run3=[];
run4=[];
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        sFrame=eset.expt(i).track(j).startFrame;
        for k=1:length(eset.expt(i).track(j).run)
             runlen=eset.expt(i).track(j).run(k).runTime;
             runstart=eset.expt(i).track(j).run(k).startInd; 
             runend=eset.expt(i).track(j).run(k).endInd;
             runmid=(runstart+runend)/2;
             runmid=fix(runmid);
             runf=sFrame+runmid;
             if runlen>50 
                 tracktheta=eset.expt(i).track(j).getDerivedQuantity('theta');
                 runtheta=tracktheta(runmid-10:runmid+10);
                 runtheta=unwrap(runtheta);
                 meanruntheta=mean(runtheta);
                 if meanruntheta>pi
                     meanruntheta=meanruntheta-2*pi;
                 elseif meanruntheta<-pi
                     meanruntheta=meanruntheta+2*pi;
                 end
                 rspeed=eset.expt(i).track(j).getDerivedQuantity('speed'); %in pixle per s
                 rspeed=rspeed/8.836;
                 meanspeed=mean(rspeed(runmid-10:runmid+10));
                 if meanruntheta>pi*3/4 || meanruntheta<-pi*3/4
                     run1=[run1, meanspeed];
                 elseif meanruntheta>pi/4 && meanruntheta<pi*3/4
                     run2=[run2, meanspeed];
                 elseif meanruntheta>-pi/4 && meanruntheta<pi/4
                     run3=[run3, meanspeed];
                 elseif meanruntheta>-pi*3/4 && meanruntheta<-pi/4
                     run4=[run4, meanspeed];
                 end
             end

        end
    end
end

m1=mean(run1);
e1=std(run1)/sqrt(length(run1));
m2=mean(run2);
e2=std(run2)/sqrt(length(run2));
m3=mean(run3);
e3=std(run3)/sqrt(length(run3));
m4=mean(run4);
e4=std(run4)/sqrt(length(run4));

speed=[m3  m1 ];
espeed=[e3  e1 ];
figure;
bar(speed);
hold on;
errorbar(speed,espeed,'r*');

                