function [thetaMatrix3D,thetaHeaders,thetaFileList,thetaDist3D,thetaDistHeaders] = JoshGetTheta_yaml(eset)


thetaMatrix3D =[]; %combined list of runs
    % need to deal with issue of different number of runs for different
    % tracks & experiments
thetaFileList = {}; % input file for each run
thetaHeaders = {};
    
thetaDist3D = []; %matrix with values for attributes binned by angle of run (theta binned)
                % Composed of: pathLength array, Xdisp array, frequency
                % array (binCount)
thetaDistHeaders = {};
binCount = []; % binCount is array of theta binned frequencies

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
TotalSecondRev = 0;
curveIndex = 0;
runLength = 0;
velocity = 0;




for i = 1:20
    bins = [bins (bin_idx*i-(bin_idx/2))];
end

meanThetas = [];

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
   
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
                pathLength = eset(i).expt(j).track(1,k).run(1,l).pathLength;
                
               % if k == 14 && l == 16
                %display(k);    
               % display(l);
               % end    
                


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
                % fprintf(textfile,'%s,%i,%i,%f,%f,%f,%f,%i,%f,%f,%f \n',eset(i).expt(j).fname,k,l,pathLength,meanTheta,startTheta,endTheta,isotherm,runTime,curveIndex,velocity);
                
                % compile run paramater values into expt matrix, thetaMatrix
                thetaFileList{l,k,j} = eset(i).expt(j).fname; % row(l)=run number, column(k)=track, z(j)=experiment
                thetaMatrix3D(end+1,:,j) = [k,l,pathLength,meanTheta,startTheta,endTheta,isotherm,runTime,curveIndex,velocity];
                %need to deal with distinction between zero and empty values better... 
                thetaHeaders = {'k','l','pathLength','meanTheta','startTheta','endTheta','isotherm','runTime','curveIndex','velocity'};
                
                % Compile frequency distribution data into matrix, thetaDist3D 
                
                [binCount,xout]=hist(meanThetas, bins);
                binCountY=transpose(binCount);
                xoutY=transpose(xout);
                paths(ceil(abs(meanTheta)/bin_idx),nnz(paths(ceil(abs(meanTheta)/bin_idx),:))+1)=pathLength;
                disps(ceil(abs(meanTheta)/bin_idx),nnz(disps(ceil(abs(meanTheta)/bin_idx),:))+1)=abs(eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).endInd)-eset(i).expt(j).track(1,k).dq.iloc(1,eset(i).expt(j).track(1,k).run(1,l).startInd));
        
                dispAvg=zeros(20,1);
                for m = 1:20
                    dispAvg(m) = sum(disps(m,:))/nnz(disps(m,:));
                end
                
                pathLAvg=zeros(20,1);
                for m = 1:20
                    pathLAvg(m) = sum(paths(m,:))/nnz(paths(m,:));
                end
                
                % name of files
                thetaDistHeaders = {'bins','Count','displacement','path length'};
           
                thetaDist3D(:,1,j)=[xoutY];
                thetaDist3D(:,2,j)=[binCountY];                 
                thetaDist3D(:,3,j)=[dispAvg];
                thetaDist3D(:,4,j)=[pathLAvg];                
                
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