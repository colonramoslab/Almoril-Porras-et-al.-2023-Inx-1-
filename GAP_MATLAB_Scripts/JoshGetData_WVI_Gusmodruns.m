numSharpTurns = 0;
totalRunLength = 0;
averageRunLength = 0;
finalXPosition = 0;
startingXPosition = 0;
posRunLength = 0;
XDisp = 0;
posRuns = 0;
posBiasCount = 0;  % will use to create a frequency bias weighted by direction vector
negBiasCount = 0;  % will use to create a frequency bias weighted by direction vector
negRuns = 0;
negXDisp = 0;
posXDisp = 0;
negRunLength = 0;
ttxi = 0;
ttxv = 0;
xdispi=0;
turnOmega=0;
turnBlip=0;
turnReversal=0;
turnSecondRev=0;
l=0;
FracNeg = 0;
NegLengthBias = 0;
posRunTime = 0;
posVelocity = 0;
posCurveIndex = 0;
negRunTime = 0;
negVelocity = 0;
negCurveIndex = 0;
WVIndex = 0;
posWVIndex = 0;
negWVIndex = 0;
TrackDuration = 0;

posStartRunLength = 0;
negStartRunLength = 0;
posStartRuns = 0;
negStartRuns = 0;
posStartBiasCount = 0;
negStartBiasCount = 0;

%outputFile = input('Type address of output file:', 's');

if exist('textfile','var')
else
    textfile=1;
end

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
        % fprintf(textfile,'%s\n',eset(i).expt(j).fname); MOVED to each
        % line to improve compiling multiple files with uniquely
        % identifying info on each line
        fprintf(textfile,'Expt,Track,Run, Run starting X Position, Run Duration, Run Y displacement, Run X displacement, Run Total displacement, Trayectory Length \n');
            %%Finding the max and min values of all tracks
            fignum=figure();
            clf(fignum);  %create a figure and clear it
            hold on;
            axis([0 2500 0 2000]);
            yL = get(gca,'YLim');
            minXPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(1);
            minYPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(2);
            maxXPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(1);
            maxYPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(2);
            for k = 1:length(eset(i).expt(j).track)
                for l = 1:length(eset(i).expt(j).track(1,k).pt)
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(1)<minXPosition
                        minXPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(1);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(2)<minYPosition
                        minYPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(2);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(1)>maxXPosition
                        maxXPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(1);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(2)>maxYPosition
                        maxYPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(2);
                    end
                end
            end
            areacheck = (maxXPosition-minXPosition)*(maxYPosition-minYPosition);
            if areacheck > 828100? %910px x 910px area
                    warning(['Tracks found in an area larger than estimated arena']);
            elseif areacheck < 720801 %849px x 849px area
                    warning(['Tracks found in an area much smaller than estimated arena']);
            end
            

        for k = 1:length(eset(i).expt(j).track);
            
            %Get distance traveled (x displacement)
            startingXPosition = eset(i).expt(j).track(1,k).pt(1,1).imOffset(1);
            finalXPosition = eset(i).expt(j).track(1,k).pt(1,length(eset(i).expt(j).track(1,k).pt)).imOffset(1);
            trackLength = 0;
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                trackLength = trackLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
            end
                if trackLength<100 %arbitrary
                    plotPath(eset(i).expt(j).track(k),'sloc','k-')
                end
                if startingXPosition < (minXPosition+10) %10px from minX
                    plotPath(eset(i).expt(j).track(k))
                end
                if startingXPosition > (maxXPosition-10) %10px from maxX
                    plotPath(eset(i).expt(j).track(k),'sloc','r-')
                end


            %Calculate average run length for all runs & each run direction
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                runStart = eset(i).expt(j).track(1,k).run(1,l).startInd;
                runEnd = eset(i).expt(j).track(1,k).run(1,l).endInd;
                totalRunLength = totalRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                % runDuration = (runEnd - runStart)/2; %% Duration of the run in s, does not take into account acquisition is not a perfect 2fps
                runTime = eset(i).expt(j).track(1,k).run(1,l).runTime; %% Pulled from JoshGetTheta
                
                %Calculations for positive (meanTheta<pi/2) & negative (meanTheta>pi/2) runs performed separately
                if abs(eset(i).expt(j).track(1,k).run(1,l).meanTheta) < (pi/2)
                   % if eset(i).expt(j).track(1,k).run(1,l).endInd > eset(i).expt(j).track(1,k).npts
                   %     display(k);
                   %     display(l);
                   % else
                        posXDisp = posXDisp + (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart));
                   % end
  
                    % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                    X0Disp = (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart)); %individual run X displacement
                    Y0Disp = (eset(i).expt(j).track(1,k).dq.iloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(2, runStart)); %individual run Y displacement
                    T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                    startingX0Position = eset(i).expt(j).track(1,k).dq.iloc(1,runStart); %individual run X starting position
                    fprintf(textfile,'%s,%i,%i,%f,%f,%f,%f,%f,%f \n',eset(i).expt(j).fname,k,l,startingX0Position, runTime, Y0Disp, X0Disp, T0Disp, totalRunLength);
                else
                   % if eset(i).expt(j).track(1,k).run(1,l).endInd > eset(i).expt(j).track(1,k).npts
                   %     display(k);
                   %     display(l);
                   % else
                        negXDisp = negXDisp + (eset(i).expt(j).track(1,k).dq.iloc(1,runStart)-eset(i).expt(j).track(1,k).dq.iloc(1,runEnd));
                   % end

                    % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                    X0Disp = (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart)); %individual run X displacement
                    Y0Disp = (eset(i).expt(j).track(1,k).dq.iloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(2,runStart)); %individual run Y displacement
                    T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                    startingX0Position = eset(i).expt(j).track(1,k).dq.iloc(1,runStart); %individual run X starting position
                    fprintf(textfile,'%s,%i,%i,%f,%f,%f,%f,%f,%f \n',eset(i).expt(j).fname,k,l,startingX0Position, runTime, Y0Disp, X0Disp, T0Disp, totalRunLength);
            % added "posTimeAvg, negTimeAvg, posVelocityAvg, negVelocityAvg, posCurveIndexAvg, negCurveIndexAvg, posWVIAvg, negWVIAvg,StartCountBias,StartLengthBias" before ttxi
            % also included expt name in each line
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                
            end
            totalRunLength = 0;
            posXDisp = 0;
            negXDisp = 0;
            posRunLength = 0;
            negRunLength = 0;
            posRuns = 0;
            negRuns = 0;
            turnOmega=0;
            turnBlip=0;
            turnReversal=0;
            turnSecondRev=0;
            posBiasCount = 0;   
            negBiasCount = 0; 
            posRunTime = 0;            
            negRunTime = 0;
            posVelocity = 0;
            negVelocity = 0;
            posCurveIndex = 0;
            negCurveIndex = 0;
            WVIndex = 0;
            posWVIndex = 0;
            negWVIndex = 0;
                end
            end
        end
    end
end
