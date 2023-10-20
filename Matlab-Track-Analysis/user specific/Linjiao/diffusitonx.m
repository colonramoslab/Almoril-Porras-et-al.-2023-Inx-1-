%calculate drift velocity/speed for salt plasticity 
xspeed=[];
t=[];
% for i=1:length(eset.expt) % all time
%     for j=1:length(eset.expt(i).track)
%         %if eset.expt(i).track(j).npts>600 && eset.expt(i).track(j).startFrame<120
%            location=eset.expt(i).track(j).getDerivedQuantity('iloc');
%            xk=location(1,:);
%            tk=eset.expt(i).track(j).getDerivedQuantity('eti');
%            xbb=xk-xk(1);
%            xbb=xbb/8.836; % to mm per sec
%            tbb=tk-tk(1);
%            %xend=min(length(xbb),1200);
%            %xp=xbb(xend)/tbb(xend);
%            xp=xbb(end)/tbb(end);
%            xspeed=[xspeed, xp];
%            t=[t,tbb(end)];
%         %end
%     end
% end


for i=1:length(eset.expt) % first 15 min
    for j=1:length(eset.expt(i).track)
        if eset.expt(i).track(j).startFrame<1800
           location=eset.expt(i).track(j).getDerivedQuantity('iloc');
           xk=location(1,:);
           tk=eset.expt(i).track(j).getDerivedQuantity('eti');
           xbb=xk-xk(1);
           xbb=xbb/8.836; % to mm per sec
           tbb=tk-tk(1);
           xend=min(length(xbb),1800);
           xp=xbb(xend)/tbb(xend);
           xp=xbb(end)/tbb(end);
           xspeed=[xspeed, xp];
           t=[t,tbb(end)];
        end
    end
end

runsp = eset.gatherFromSubField('run', 'speed');
runsp=runsp/8.836;
wmean(xspeed,t)
sqrt(var(xspeed,t)/(20*i))
mean(runsp)
std(runsp)/sqrt(i)
xspeed=xspeed/mean(runsp);
wmean(xspeed,t)
sqrt(var(xspeed,t)/(20*i))