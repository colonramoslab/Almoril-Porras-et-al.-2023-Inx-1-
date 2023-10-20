num=100;

for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        location=eset.expt(i).track(j).getDerivedQuantity('sloc');
        sFrame=eset.expt(i).track(j).startFrame;
        for k=1:length(eset.expt(i).track(j).reorientation)
            reo=eset.expt(i).track(j).reorientation(k);
            reoseq=reo.turnsequence;
%            if isequal(reoseq, -1)   %omega turn
            if isequal(reoseq, [1 -1])   %reversal omega turn
                run1=eset.expt(i).track(j).run(k);
                run2=eset.expt(i).track(j).run(k+1);
                in=run1.endTheta;
                out=run2.startTheta;
                ind1=run1.endInd;
                ind2=run2.startInd;
                ind0=reo.sharpTurn(1).centralInd;
                delta=reo.dTheta;
                delta=rad2deg(delta);
                x01=location(1,ind1);
                y01=location(2,ind1);
                x02=location(1,ind2);
                y02=location(2,ind2);
                x0=location(1,ind0);
                y0=location(2,ind0);
%                 x1=location(1,ind1-5:ind1)-x01;
%                 y1=location(2,ind1-5:ind1)-y01;
%                 x2=location(1,ind2:ind2+5)-x02;
%                 y2=location(2,ind2:ind2+5)-y02;
                x1=location(1,ind1-5:ind0)-x0;
                y1=location(2,ind1-5:ind0)-y0;
                x2=location(1,ind0:ind2+5)-x0;
                y2=location(2,ind0:ind2+5)-y0;
                newx1=x1*cos(in)+y1*sin(in);
                newy1=-x1*sin(in)+y1*cos(in);
                newx2=x2*cos(in)+y2*sin(in);
                newy2=-x2*sin(in)+y2*cos(in);
                %figure
                if num>0 && (sFrame+ind1)<1800
                    num=num-1;
%                  plot(newx1,newy1);
                hold on;
                  plot(newx2,newy2);
                if isequal(reoseq, [1 -1])
                    indO=reo.sharpTurn(2).centralInd;
                    xO=location(1,indO)-x0;
                    yO=location(2,indO)-y0;
                    newxO=xO*cos(in)+yO*sin(in);
                    newyO=-xO*sin(in)+yO*cos(in);
                    plot(newxO,newyO,'ms','MarkerSize',5)
                end
                    if num==2
                      plot(newx2,newy2,'g');
                    plot(newxO,newyO,'gs','MarkerSize',5)
%                       plot(newx1,newy1,'g');
                    end            
                end
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

hold on;
plot(0,0,'m.','MarkerSize',20);
axis equal;

                