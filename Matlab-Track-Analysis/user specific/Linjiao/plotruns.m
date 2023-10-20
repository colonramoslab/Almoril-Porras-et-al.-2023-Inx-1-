num=100;

for i=3:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        location=eset.expt(i).track(j).getDerivedQuantity('sloc');
        sFrame=eset.expt(i).track(j).startFrame;
        for k=1:length(eset.expt(i).track(j).reorientation)
             runlen=eset.expt(i).track(j).run(k).runTime;
             runstart=eset.expt(i).track(j).run(k).startInd; 
             runend=eset.expt(i).track(j).run(k).endInd;
             runmid=(runstart+runend)/2;
             runmid=fix(runmid);
             runf=sFrame+runmid;
             if runlen>50 && runf<1800 && num>0
                 num=num-1;
                 tracktheta=eset.expt(i).track(j).getDerivedQuantity('theta');
                 runtheta=tracktheta(runmid-10:runmid+10);
                 d_run=runtheta(end)-runtheta(1);
                 runtheta=unwrap(runtheta);
                 %meanruntheta=mean(runtheta);
                 meanruntheta=mean(runtheta(1:5));
                 x0=location(1,runmid-10);
                 y0=location(2,runmid-10);
                 x=location(1,(runmid-10):(runmid+10))-x0;
                 y=location(2,(runmid-10):(runmid+10))-y0;
                 newx=x*cos(meanruntheta)+y*sin(meanruntheta);
                 newy=-x*sin(meanruntheta)+y*cos(meanruntheta);
                 if abs(d_run)<pi/2
                 plot(newx,newy);
                 hold on;
                 end
             end

        end
    end
end

hold on;
plot(0,0,'m','MarkerSize',20);
axis equal;

                