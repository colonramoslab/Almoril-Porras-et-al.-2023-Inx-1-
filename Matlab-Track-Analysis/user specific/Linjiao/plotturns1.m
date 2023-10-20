
for i=1:length(eset1.expt)
    for j=1:length(eset1.expt(i).track)
        location=eset1.expt(i).track(j).getDerivedQuantity('sloc');
        for k=1:length(eset1.expt(i).track(j).reorientation)
            reo=eset1.expt(i).track(j).reorientation(k);
            reoseq=reo.turnsequence;
%            if isequal(reoseq, -1)   %omega turn
            if isequal(reoseq, [1 -1])   %reversal omega turn
                run1=eset1.expt(i).track(j).run(k);
                run2=eset1.expt(i).track(j).run(k+1);
                in=run1.endTheta;
                out=run2.startTheta;
                ind1=run1.endInd;
                ind2=run2.startInd;
                delta=reo.dTheta;
                delta=rad2deg(delta);
                x01=location(1,ind1);
                y01=location(2,ind1);
                x02=location(1,ind2);
                y02=location(2,ind2);
                x1=location(1,ind1-10:ind1)-x01;
                y1=location(2,ind1-10:ind1)-y01;
                x2=location(1,ind2:ind2+10)-x02;
                y2=location(2,ind2:ind2+10)-y02;
                newx1=x1*sin(in)-y1*cos(in);
                newy1=x1*cos(in)+y1*sin(in);
                newx2=x2*sin(in)-y2*cos(in);
                newy2=x2*cos(in)+y2*sin(in);
                %figure
                %plot(newx1,newy1,'r.-');
                hold on;
                plot(0,0,'m.','MarkerSize',20);
                plot(newx2,newy2,'.-');
                axis equal;
%                 title (['number of reorientation =' num2str(i) ' ' num2str(j) ' ' num2str(k) ' dtheta=' num2str(delta)]);
%                 figure;
%                 eset.expt(i).track(j).plotPath('sloc','.-','inds',s:e);
%                 hold on;
%                 eset.expt(i).track(j).plotPath('sloc','m*','inds',c:c);
%                 axis equal;
            end
        end
    end
end
                