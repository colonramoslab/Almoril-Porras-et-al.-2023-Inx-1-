pathLength = 0;
meanTheta = 0;
startTheta = 0;
endTheta = 0;
l=0;
isotherm = 0;
runTime = 0;
bin_idx = (maxXPosition-minXPosition)/20;
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
runBin = 0;

if exist('textfile','var')
else
    textfile=1;
end

for i = 1:20
    bins = [bins (bin_idx*i-(bin_idx/2))];
end

binrange = [];
for i = 1:21
    binrange = [binrange (minXPosition+bin_idx*(i-1))];
end

meanThetas = [];

for i = 1:length(eset)
    for j = 1:length(eset(i).expt)
        % fprintf(textfile,'%s\n',eset(i).expt(j).fname); % Moved to each
        % line on 3/23.13
        fprintf(textfile,'Expt,Track,Run,Isotherm,Run Length,Run Time,Start X Position,curveIndex,velocity,Mean Theta,Start Theta,End Theta\n');
        for k = 1:length(eset(i).expt(j).track)
            
            %Calculate average run length
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                meanTheta = eset(i).expt(j).track(1,k).run(1,l).meanTheta;
                if meanTheta == 0
                    meanTheta = 0.0001;
                end  
                meanThetas = [meanThetas abs(meanTheta)];
                startTheta = eset(i).expt(j).track(1,k).run(1,l).startTheta;
                endTheta = eset(i).expt(j).track(1,k).run(1,l).endTheta;
                pathLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                paths(ceil(abs(meanTheta)/bin_idx),nnz(paths(ceil(abs(meanTheta)/bin_idx),:))+1)=pathLength;
                disps(ceil(abs(meanTheta)/bin_idx),nnz(disps(ceil(abs(meanTheta)/bin_idx),:))+1)=abs(eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).endInd)-eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).startInd));
                runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                
                                                        % !!! NEW CODE BY JDH, 03/22/13 !!! %
                % Define run indices in coordinate array, iloc                                  
                runStart = eset(i).expt(j).track(1,k).run(1,l).startInd;
                runEnd = eset(i).expt(j).track(1,k).run(1,l).endInd;
                % calculate run displacement, sqrt(xdisp^2+ydisp^2)
                X0Disp = (eset(i).expt(j).track(1,k).dq.sloc(1,runEnd)-eset(i).expt(j).track(1,k).dq.sloc(1,runStart)); %individual run X displacement
                Y0Disp = (eset(i).expt(j).track(1,k).dq.sloc(2,runEnd)-eset(i).expt(j).track(1,k).dq.sloc(2,runStart)); %individual run Y displacement
                T0Disp = sqrt(X0Disp^2 + Y0Disp^2); %individual run total displacement vector length
                % calculate curvature index, Displacement/Run length,
                runLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                curveIndex = T0Disp/runLength;
                % calculate run velocity, displacement/time
                % runTime = eset(i).expt(j).track(1,k).run(1,l).runTime;
                velocity = T0Disp/runTime;
                
                %New Code by GAP 2020/08/06
                StartX0 = eset(i).expt(j).track(1,k).dq.sloc(1,runStart);
                
                if abs(X0Disp)<11.25
                    if abs(Y0Disp/T0Disp)>0.9
                    isotherm = 1;
                    end
                end
                
                % added filename at beginning of each line on 3/23/13
                % added curve index & velocity to output on 3/23
                fprintf(textfile,'%s,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f \n',eset(i).expt(j).fname,k,l,isotherm,pathLength,runTime,StartX0,curveIndex,velocity,meanTheta,startTheta,endTheta);
                
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
    end
end
