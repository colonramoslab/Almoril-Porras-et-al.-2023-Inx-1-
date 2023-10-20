pathLength = 0;
meanTheta = 0;
startTheta = 0;
endTheta = 0;
l=0;
isotherm = 0;
runTime = 0;
bin_idx = pi/20;
bins = [];
paths = zeros(20,0);
sparse(paths);
disps = zeros(20,0);
sparse(disps);
turnOmega = 0;
turnBlip = 0;
turnReversal = 0;
turnSecondRev = 0;
%outputFile = input('Type address of output file:', 's');
TotalOmega = 0;
TotalBlip = 0;
TotalReversal = 0;
TotalSecondRev = 0;
curveIndex = 0;
runLength = 0;
velocity = 0;

if exist('textfile','var')
else
    textfile=1;
end

for i = 1:20
    bins = [bins (bin_idx*i-(bin_idx/2))];
end

meanThetas = [];

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
        % fprintf(textfile,'%s\n',eset(i).expt(j).fname); % Moved to each
        % line on 3/23.13
        fprintf(textfile,'Expt, Track,Run,Run Length,Mean Theta,Start Theta,End Theta,Isotherm,Run Time,curveIndex,velocity\n');
        for k = 1:length(eset(i).expt(j).track);            
            
            %Calculate average run length
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                meanTheta = eset(i).expt(j).track(1,k).run(1,l).meanTheta;
                if meanTheta == 0
                    meanTheta = 0.0001;
                end  
                meanThetas = [meanThetas abs(meanTheta)];
                startTheta = eset(i).expt(j).track(1,k).run(1,l).startTheta;
                endTheta = eset(i).expt(j).track(1,k).run(1,l).endTheta;
                
               % if k == 14 && l == 16
                %display(k);    
               % display(l);
               % end    
                
                pathLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                paths(ceil(abs(meanTheta)/bin_idx),nnz(paths(ceil(abs(meanTheta)/bin_idx),:))+1)=pathLength;
                disps(ceil(abs(meanTheta)/bin_idx),nnz(disps(ceil(abs(meanTheta)/bin_idx),:))+1)=abs(eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).endInd)-eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).startInd));
                runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                if abs(cos(meanTheta))< 0.14
                    isotherm = 1;
                end
                
                                                        % !!! NEW CODE BY JDH, 03/22/13 !!! %
                % Define run indices in coordinate array, iloc                                  
                runStart = eset(i).expt(j).track(1,k).run(1,l).startInd;
                runEnd = eset(i).expt(j).track(1,k).run(1,l).endInd;
                % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                X0Disp = (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart)); %individual run X displacement
                Y0Disp = (eset(i).expt(j).track(1,k).dq.iloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(2,runStart)); %individual run Y displacement
                T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                % calculate curvature index, Displacement/Run length,
                runLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                curveIndex = T0Disp/runLength;
                % calculate run velocity, displacement/time
                % runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                velocity = T0Disp/runTime;
                
                % added filename at beginning of each line on 3/23/13
                % added curve index & velocity to output on 3/23
                fprintf(textfile,'%s,%i,%i,%f,%f,%f,%f,%i,%f,%f,%f \n',eset(i).expt(j).fname,k,l,pathLength,meanTheta,startTheta,endTheta,isotherm,runTime,curveIndex,velocity);
                
            pathLength = 0;
            meanTheta = 0;
            startTheta = 0;
            endTheta = 0;
            isotherm = 0;
            runTime = 0;
            curveIndex= 0;
            runLength = 0;
            velocity = 0;
            end
            
            
            
            
        end
       
        
        %Print histogram of average theta values
        [n,xout]=hist(meanThetas, bins);
        xout = xout * -1;
        fignum=figure();
        clf(fignum);  %create a figure and clear it
        hold on;
        bar(xout(1:10),n(1:10),'r');
        bar(xout(11:20),n(11:20),'b');
        if isnan(max(n))||max(n)==0
            ul=1;
        else
            ul=ceil(max(n)*1.2);
        end
        axis([-pi 0 0 ul]);
        yL = get(gca,'YLim');
        line([-pi/2 -pi/2],yL,'Color','k');
        hold off;
        set(gcf, 'Position', get(0,'Screensize'));
        saveas(fignum,strcat(outputdirectory,fileprefix,'_theta_hist.jpg'));
        close(fignum);
        
        %Print average path length grouped by theta values
        %foo1=mean(paths,2);
        
        foo1=zeros(20,1);
        for m = 1:20
            foo1(m) = sum(paths(m,:))/nnz(paths(m,:));
        end
        fignum=figure();
        clf(fignum);  %create a figure and clear it
        hold on;
        bar(xout(1:10),foo1(1:10),'r');
        bar(xout(11:20),foo1(11:20),'b');
        %bar(xout,foo1);
        if isnan(max(foo1))||max(foo1)==0
            ul=1;
        else
            ul=ceil(max(foo1)*1.2);
        end
        axis([-pi 0 0 ul]);
        yL = get(gca,'YLim');
        line([-pi/2 -pi/2],yL,'Color','k');
        hold off;
        set(gcf, 'Position', get(0,'Screensize'));
        saveas(fignum,strcat(outputdirectory,fileprefix,'_theta_path_avg.jpg'));
        close(fignum);
        
        %Print average path length grouped by theta values
        %foo=mean(disps,2);
        foo=zeros(20,1);
        for m = 1:20
            foo(m) = sum(disps(m,:))/nnz(disps(m,:));
        end
        fignum=figure();
        clf(fignum);  %create a figure and clear it
        hold on;
        bar(xout(1:10),foo(1:10),'r');
        bar(xout(11:20),foo(11:20),'b');
        %bar(xout,foo);
        if isnan(max(foo1))||max(foo1)==0
            ul=1;
        else
            ul=ceil(max(foo1)*1.2);
        end
        axis([-pi 0 0 ul]);
        yL = get(gca,'YLim');
        line([-pi/2 -pi/2],yL,'Color','k');
        hold off;
        set(gcf, 'Position', get(0,'Screensize'));
        saveas(fignum,strcat(outputdirectory,fileprefix,'_theta_XDisp_avg.jpg'));
        close(fignum);
        
        fprintf(textfile,'\n\nBin,Path Length Average,X Disp Average,Count\n');
        
        for v = 1:20
            fprintf(textfile,'%f,%f,%f,%i\n',xout(v),foo1(v),foo(v),n(v));
        end
        
        
    end
end
