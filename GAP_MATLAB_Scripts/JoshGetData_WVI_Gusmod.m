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
        fprintf(textfile,'Expt,Track,Run,# Sharp Turns,Omega Turns,Blip Turns,Reversals,Second Reversals,Distance Traveled X,Average Run Length,Positive Average Run Length,Positive Average X Displacement,Negative Average Run Length,Negative Average X Displacement,Positive Run Count,Negative Run Count,Fraction Neg Runs,Negative Bias Count,Positive Bias Count,posTimeAvg,negTimeAvg,posVelocityAvg,negVelocityAvg,posCurveIndexAvg,negCurveIndexAvg,posWVIAvg,negWVIAvg,StartCountBias,StartLengthBias,TTXI,TTXV,absTTXV,starting X Position,Track Duration \n');

        for k = 1:length(eset(i).expt(j).track);
            %Get # of sharp turns
            numSharpTurns = length(eset(i).expt(j).track(1,k).sharpTurn);
            
            for m = 1:length(eset(i).expt(j).track(1,k).sharpTurn)
                if eset(i).expt(j).track(1,k).sharpTurn(1,m).typeCode == -1
                    turnOmega = turnOmega + 1;
                elseif eset(i).expt(j).track(1,k).sharpTurn(1,m).typeCode == 0
                    turnBlip = turnBlip + 1;
                elseif eset(i).expt(j).track(1,k).sharpTurn(1,m).typeCode == 1
                    turnReversal = turnReversal + 1;
                elseif eset(i).expt(j).track(1,k).sharpTurn(1,m).typeCode == 2
                    turnSecondRev = turnSecondRev +1;
                end
            end
            
            %Get distance traveled (x displacement)
            startingXPosition = eset(i).expt(j).track(1,k).pt(1,1).imOffset(1);
            finalXPosition = eset(i).expt(j).track(1,k).pt(1,length(eset(i).expt(j).track(1,k).pt)).imOffset(1);
            TrackDuration = (eset(i).expt(j).track(1,k).endFrame - eset(i).expt(j).track(1,k).startFrame)/2; %Duration of the track in s
            
            %Calculate average run length for all runs & each run direction
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                runStart = eset(i).expt(j).track(1,k).run(1,l).startInd;
                runEnd = eset(i).expt(j).track(1,k).run(1,l).endInd;
                totalRunLength = totalRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                
                %Calculations for positive (meanTheta<pi/2) & negative (meanTheta>pi/2) runs performed separately
                if abs(eset(i).expt(j).track(1,k).run(1,l).meanTheta) < (pi/2)
                   % if eset(i).expt(j).track(1,k).run(1,l).endInd > eset(i).expt(j).track(1,k).npts
                   %     display(k);
                   %     display(l);
                   % else
                        posXDisp = posXDisp + (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart));
                   % end
                    posRunLength = posRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                    posRuns = posRuns + 1;
                    posBiasCount = posBiasCount + cos(eset(i).expt(j).track(1,k).run(1,l).meanTheta);
  
                    % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                    X0Disp = (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart)); %individual run X displacement
                    Y0Disp = (eset(i).expt(j).track(1,k).dq.iloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(2, runStart)); %individual run Y displacement
                    T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                    % calculate curvature index, Displacement/Run length,
                    runLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                    curveIndex = T0Disp/runLength;
                    % calculate run velocity, displacement/time
                    runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                    velocity = T0Disp/runTime;
                    % accumulate measures for average
                    posRunTime = posRunTime + runTime;
                    posVelocity = posVelocity + velocity;
                    posCurveIndex = posCurveIndex + curveIndex;
                    % calculate weathervane index, WVI = 
                    % (abs(end theta)-abs(start theta))/pi+1)/2   
                    % For cold-directed weathervaning, 0.5<WVI<1.0
                    % For warm-directed weathervaining, 0<WVI<0.5
                    startTheta = eset(i).expt(j).track(1,k).run(1,l).startTheta;
                    endTheta = eset(i).expt(j).track(1,k).run(1,l).endTheta;
                    WVIndex = ((abs(endTheta)-abs(startTheta))/pi+1)/2;
                    posWVIndex = posWVIndex + WVIndex;
                   
                    
                else
                   % if eset(i).expt(j).track(1,k).run(1,l).endInd > eset(i).expt(j).track(1,k).npts
                   %     display(k);
                   %     display(l);
                   % else
                        negXDisp = negXDisp + (eset(i).expt(j).track(1,k).dq.iloc(1,runStart)-eset(i).expt(j).track(1,k).dq.iloc(1,runEnd));
                   % end
                    negRunLength = negRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                    negRuns = negRuns + 1;
                    negBiasCount = negBiasCount + cos(eset(i).expt(j).track(1,k).run(1,l).meanTheta);
                    
                    % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                    X0Disp = (eset(i).expt(j).track(1,k).dq.iloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(1,runStart)); %individual run X displacement
                    Y0Disp = (eset(i).expt(j).track(1,k).dq.iloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.iloc(2,runStart)); %individual run Y displacement
                    T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                    % calculate curvature index, Displacement/Run length,
                    runLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                    curveIndex = T0Disp/runLength;
                    % calculate run velocity, displacement/time
                    runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                    velocity = T0Disp/runTime;
                    % accumulate measures for average
                    negRunTime = negRunTime + runTime;
                    negVelocity = negVelocity + velocity;
                    negCurveIndex = negCurveIndex + curveIndex;
                    % calculate weathervane index of course correction  
                    %           WVI= (abs(end theta)-abs(start theta))/pi+1)/2   
                    % For cold-directed weather-vaning, 0.5<WVI<1.0
                    % For warm-directed weather-vaning, 0<WVI<0.5
                    startTheta = eset(i).expt(j).track(1,k).run(1,l).startTheta;
                    endTheta = eset(i).expt(j).track(1,k).run(1,l).endTheta;
                    WVIndex = ((abs(endTheta)-abs(startTheta))/pi+1)/2;
                    negWVIndex = negWVIndex + WVIndex;
                end
                
                % collect statistics on runs that *begin* in warm direction    
                if abs(eset(i).expt(j).track(1,k).run(1,l).startTheta) < (pi/2)
                        posStartRunLength = posStartRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                        posStartRuns = posStartRuns + 1;
                        posStartBiasCount = posStartBiasCount + cos(eset(i).expt(j).track(1,k).run(1,l).meanTheta);
                else
                        negStartRunLength = negStartRunLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
                        negStartRuns = negStartRuns + 1;
                        negStartBiasCount = negStartBiasCount + cos(eset(i).expt(j).track(1,k).run(1,l).meanTheta);
                end
                
             end
                       
                   
            
            averageRunLength = totalRunLength / l;
           
            if posRuns == 0
                   posRunAvg = 0;
                   posXDispAvg = 0;
                   posTimeAvg = 0;  % JDH, 3/22/13
                   posVelocityAvg = 0;  % JDH, 3/22/13
                   posCurveIndexAvg = 0;  % JDH, 3/22/13
                   posWVIAvg = 0;
            else
                posRunAvg = posRunLength/posRuns;
                posXDispAvg = (cast(posXDisp,'double')/posRuns);
                posTimeAvg = posRunTime/posRuns;  % JDH, 3/22/13
                posVelocityAvg = posVelocity/posRuns;  % JDH, 3/22/13
                posCurveIndexAvg = posCurveIndex/posRuns; % JDH, 3/22/13
                posWVIAvg=posWVIndex/posRuns;
            end
            
             if negRuns == 0
                   negRunAvg = 0;
                   negXDispAvg = 0;
                   negTimeAvg = 0;  % JDH, 3/22/13
                   negVelocityAvg = 0;  % JDH, 3/22/13
                   negCurveIndexAvg = 0;  % JDH, 3/22/13
                   negWVIAvg = 0;
             else
                negRunAvg = negRunLength/negRuns;
                negXDispAvg = (cast(negXDisp,'double')/negRuns);
                negTimeAvg = negRunTime/negRuns;  % JDH, 3/22/13
                negVelocityAvg = negVelocity/negRuns;  % JDH, 3/22/13
                negCurveIndexAvg = negCurveIndex/negRuns; % JDH, 3/22/13
                negWVIAvg = negWVIndex/negRuns;
             end
             
            if posStartRuns == 0
                   posStartRunAvg = 0;
            else
                posStartRunAvg = posStartRunLength/posStartRuns;
            end
            
             if negStartRuns == 0
                   negStartRunAvg = 0;
             else
                negStartRunAvg = negStartRunLength/negStartRuns;
             end
             
            StartCountBias = negStartRuns/(posStartRuns+negStartRuns);
            StartLengthBias = negStartRunAvg/(posStartRunAvg+negStartRunAvg);

            FracNeg = negRuns/(posRuns+negRuns);
            ttxi = (posRunAvg - negRunAvg) / (posRunAvg + negRunAvg);
            xdispi = (posXDispAvg - negXDispAvg) / (posXDispAvg + negXDispAvg);
            
            %Print output
            XDisp = (finalXPosition - startingXPosition);
            
            ttxv = (idivide((XDisp + 900),200))-4;
            if ttxv > 4
                ttxv = 4;
            elseif ttxv < -4
                ttxv = -4;
            end
            fprintf(textfile,'%s,%i,%i,%i,%i,%i,%i,%i,%f,%f,%f,%f,%f,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i,%i \n',eset(i).expt(j).fname,k,numSharpTurns,turnOmega,turnBlip,turnReversal,turnSecondRev,XDisp,averageRunLength,posRunAvg,posXDispAvg,negRunAvg,negXDispAvg,posRuns,negRuns,FracNeg,negBiasCount,posBiasCount,posTimeAvg, negTimeAvg, posVelocityAvg, negVelocityAvg, posCurveIndexAvg, negCurveIndexAvg, posWVIAvg, negWVIAvg, StartCountBias, StartLengthBias, ttxi, ttxv, abs(ttxv), startingXPosition, TrackDuration);
            % added "posTimeAvg, negTimeAvg, posVelocityAvg, negVelocityAvg, posCurveIndexAvg, negCurveIndexAvg, posWVIAvg, negWVIAvg,StartCountBias,StartLengthBias" before ttxi
            % also included expt name in each line
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
