%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MULTIPLE TEMPORAL LIGHT ANALYSIS SCRIPT
%
%   Written by Liz    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%This script compares temporal triangle and exponential ramps

% Set minimum and maximum light threshold values using target light levels
maxLT = 200;    %this should be projector target light level max that was used to aquire data
minLT = 0;      %this should be projector target light level min that was used to aquire data 
onThresh = minLT + .95*(maxLT-minLT);   %max thresh is 5% less of diff between max/min
offThresh = minLT + .05*(maxLT-minLT);  %min thresh is 5% more of diff between max/min
endSec = 1400; %cut data 23 minutes (derived by looking at all time and deciding that response was max during first 23 mins)
trimBuffer=50; %pixel buffer distance to trim within extraction window




%% LOAD EXPERIMENT
    numExpt=input('What is the number of ESETS that you want to load?');
    eset=zeros(1,numExpt);
    for j=1:length(numExpt)
        [eset(j)]=LoadFromMat;
    end
        
    if (~exist('eset','var')) || ao.resetAll  % only load if you haven't already
        matF=input('Do you want to load from MAT Files. Enter 1 or 0: ');
        if matF
            eset=LoadFromMat();
            bypass=true;
        else
            eset = ExperimentSet.fromFiles();
            bypass=false;
        end
    end
