% Script that analyzes tracks from movies of crawling larvae 
% subjected to spatial temperature gradient.
%
% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) files, is loaded into MATLAB, then processed to show how larval behavior depends on 
% temperature
% 
% START WITH SOME ANALYSIS PARAMETERS
minTrackLength = 500; % only accept tracks more than X frames in duration.
%%whichToLoad = 1; % when multiple files are in a folder, but we don't want 
                 % to look at all of them.  
lasttime = 1200; % time where we stop being interested (in seconds).
cuts = [-inf inf]; % where the offset steps are.
lengthPerPixel = 1; %pixels/sec
lengthPerPixel = 0.1054; %mm/pixel with the 8X lens; also correct boundary values
%lengthPerPixel = 0.0686; %12X lens; also correct boundary values
%  
% figureNum = 1;

% LOAD THE FILE(S)
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles();
    %use the below code to get report when analyzing data for the first
    %time
    %
    %ecl = ESetCleaner;
    %ecl.getReport(eset);
end

% CLEAN THE TRACKS
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid=.95;
    ecl.minDist=100;
    ecl.minSpeed = 1;
    ecl.minPts = minTrackLength;
    ecl.clean(eset);
    
    cleanEset = false;
    
end

% FIX HEAD-TAIL ORIENTATION, SET SEGMENT SPEEDS, DO SEGMENTATION (set
% head tail direction)
existsAndDefault('fixht',true);
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

% set segmentation threshold speeds 
existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% segment the tracks; make segment variable and do segmentation if not
% already done
existsAndDefault('segment', true);
if (segment)
    eset.executeTrackFunction('segmentTrack'); %does segmentation
    segment = false; %don't do it again when you run the script next
end


% Figure 1 (1B: thermotaxis index for each stage at each temp range)
% Set a thin zone (same one for L1/L2/L3) where x=0
positions=eset.gatherField('loc'); %gets all location pixels from all runs
xpositions=positions(1,:); %gets only x-positions
warmFrames=sum(xpositions<1297); %%For 12x lens (L1), use 1271 & 1419
coldFrames=sum(xpositions>1393); %to say "for each frame, if x><0"
thindex=(warmFrames-coldFrames)/(warmFrames+coldFrames) %use in Thindex.m file
%
% 
% Figure 2 (5K: maximum head sweep angle for each direction quadrant)
wasGoing=eset.gatherSubField('headSwing','prevDir');
% thetas=eset.gatherSubField('headSwing','prevRun','meanTheta'); %gets previous direction
maxAngle=eset.gatherSubField('headSwing','maxTheta'); %gets all maxThetas
A=maxAngle(wasGoing<(pi/4)&wasGoing>(-pi/4));
% A=A0(wasGoing>(-pi/4));
Amean=mean(A);
Ameandeg=radtodeg(Amean); 
B=maxAngle(wasGoing<(3*pi/4)&wasGoing>(pi/4));
Bmean=mean(B);
Bmeandeg=radtodeg(Bmean);
C=maxAngle(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4));
Cmean=mean(C);
Cmeandeg=radtodeg(Cmean);
D=maxAngle(wasGoing<(-pi/4)&wasGoing>(-3*pi/4));
Dmean=mean(D);
Dmeandeg=radtodeg(Dmean);
x=[1,2,3,4];   
y=[Ameandeg,Bmeandeg,Cmeandeg,Dmeandeg];
figure(5);
bar(x,y);
xlabel('direction quadrant (degrees)');
ylabel('average maximum head sweep angle (degrees)');
title('average maximum head sweep angle by direction quadrant');

% % OLD Figure 2 (5K: maximum head sweep angle for each direction quadrant)
% aveThetas=eset.gatherSubField('run','meanTheta'); %all the meanThetas for all the runs
% % MATLAB doesn't like: wasGoing=eset.gatherSubField('reorientation','prevDir'); %gets previous
% % run directions (orientation just before the start of a reorientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%when testing, matlab says error for prevDir.
%%%%%%%%%%%%%%%%%%%%%%%%%%%Checking reorientation category says
%%%%%%%%%%%%%%%%%%%%%%%%%%%prevDir=NaN. Use meanTheta instead.
% A0=aveThetas(aveThetas<(pi/4));
% A=A0(A0>(-pi/4));  %(+45 to -45 degrees) toward right of jpeg (to cold)
% Amean=mean(A);
% Ameandeg=radtodeg(Amean);
% B0=aveThetas(aveThetas<(3*pi/4));
% B=B0(B0>(pi/4)); %(+135 to +45 degrees) toward top of jpeg
% Bmean=mean(B);
% Bmeandeg=radtodeg(Bmean);
% C=aveThetas(aveThetas>(3*pi/4)|aveThetas<(-3*pi/4));
% Cmean=mean(C);
% Cmeandeg=radtodeg(Cmean);
% D0=aveThetas(aveThetas<(-pi/4));
% D=D0(D0>(-3*pi/4)); %(-135 to -45 degrees) toward bottom of jpeg
% Dmean=mean(D);
% Dmeandeg=(radtodeg(Dmean));
% x=[1,2,3,4];
% y=[Ameandeg,Bmeandeg,Cmeandeg,Dmeandeg];
% figure(1);
% bar(x,y);
% xlabel('direction quadrant (degrees)');
% ylabel('average maximum head sweep angle (degrees)');
% title('average maximum head sweep angle by direction quadrant');
%
%
% Figure 3 (5J: probability of starting new run during each head sweep when
% facing towards warmth or cold) 
accs=eset.gatherSubField('headSwing','accepted'); %gets all the "accepted" category
Aaccs=accs(wasGoing<(pi/4)&wasGoing>(-pi/4));
Anew=sum(Aaccs>0)/length(Aaccs); %(sum of all accepted headsweeps in A)/(sum of all headsweeps in A)
Caccs=accs(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4));
Cnew=sum(Caccs>0)/length(Caccs);
y=[Anew,Cnew];
figure(6); %%%%WORKS BUT NOT SURE THAT DATA MAKES SENSE (66% accept for Acool and 70% accept for Cwarm, but shouldn't more accept for Acool b/c they are trying to go away from heat and toward cool?)
bar(y);
xlabel('direction quadrant (degrees)');
ylabel('probability of starting new run during each head sweep (%)');
title('likelihood of accepting head sweep by direction quadrant');
%
%
% Figure 4 (5I: probability that head sweeps first towards warmth) 
% Input actual temperature range and midpoint. %not necessary; instead,
% have viewer interpret which is correct direction
%   if midpoint<22.5, then correct direction is toward left or x<0 (+
%   thermotaxis);
%   else correct direction is toward right or x>0 (- thermotaxis); %(adjust
%   later if it turns out different stages have different Tps).
firstmax=eset.gatherSubField('firsths','maxTheta'); %only look at
% maxThetas for the first head sweep of each reorientation
firstprevDir=eset.gatherSubField('firsths','prevDir'); %get previous run directions, but only for firstHS
warmFrames=sum((firstmax>0&firstprevDir<pi)|(firstmax<0&firstprevDir>pi)); %signs may not be right; verify
coolFrames=sum((firstmax>0&firstprevDir>pi)|(firstmax<0&firstprevDir<pi));
warmFirst=(warmFrames/(warmFrames+coolFrames))*100;
figure(7);
bar(warmFirst);
ylabel('probability head sweeps first towards warm direction (%)');
%
%
% Figure 5 (5H: number of head sweeps during each turning event)
numberHS = eset.gatherSubField('reorientation','numHS'); %numHS=number of
% HS in that reorientation
numberHS1=numberHS(numberHS>0);
figure(8);
hist(numberHS1);
% NO: numSweep=the highest numbered headswing for each reorientation where
% property valid=1; 
% numSweep1=((number of times numSweep=1)/(length(numberHS)))*100;
% numSweep2=((number of times numSweep=2)/(total number of reorientations))*100;
% numSweep3=((number of times numSweep=3)/(total number of reorientations))*100;
% x=[1,2,3];
% y=[numSweep1,numSweep2,numSweep3];
% bar(x,y);
xlabel('number of head sweeps during each turning event');
ylabel('count (%)');
%
%
% Figure 6 (5E: percent of runs in each heading direction (degrees))
% Recall from fig 2: maxAngle=eset.gatherSubField('headSwing','maxTheta'); %gets head sweep angles
A1=wasGoing(wasGoing<0&wasGoing>-pi/4); %(-45 to 0 degrees)
% toward right of jpeg (to cold)
A2=wasGoing(wasGoing<pi/4&wasGoing>0); %(0 to 45 degrees) toward
% right of jpeg (to cold)
B1=wasGoing(wasGoing<pi/2&wasGoing>pi/4); %(45 to 90 degrees) toward top of jpeg
B2=wasGoing(wasGoing<3*pi/4&wasGoing>pi/2); %(90 to 135 degrees)
% toward top of jpeg
C1=wasGoing(wasGoing<pi&wasGoing>3*pi/4); %(180 to 135 degrees) toward left of jpeg
% (to hot)
C2=wasGoing(wasGoing<(-3*pi/4)&wasGoing>-pi); %(-135 to 180 degrees) toward left of jpeg
% (to hot)
D1=wasGoing(wasGoing<(-pi/2)&wasGoing>(-3*pi/4)); %(-135 to -90 degrees) toward bottom of
% jpeg
D2=wasGoing(wasGoing<(-pi/4)&wasGoing>(-pi/2)); %(-90 to -45 degrees) toward bottom of
% jpeg
aveThetas=eset.gatherSubField('run','meanTheta'); %all the meanThetas for all the runs
numRuns=length(aveThetas); %the total number of runs
A1run=((length(A1))/numRuns)*100; 
A2run=((length(A2))/numRuns)*100;
B1run=((length(B1))/numRuns)*100;
B2run=((length(B2))/numRuns)*100;
C1run=((length(C1))/numRuns)*100;
C2run=((length(C2))/numRuns)*100;
D1run=((length(D1))/numRuns)*100;
D2run=((length(D2))/numRuns)*100;
y=[A2run,B1run,B2run,C1run,C2run,D1run,D2run,A1run];
figure(9);
bar(y); %OR could do ugly version with: hist(aveThetas,8);
xlabel('run heading (degrees)');
ylabel('count (%)');
%
%
% Figure 7 (5D: average run speed(mm/min) for each run heading)
spdraw=eset.gatherField('speed','run');
spd=spdraw*lengthPerPixel*240; %convert to mm/min
% speeds=eset.gatherSubField('run','speed','mean'); %OR do this to get mean
% SPEEDS for each run %this line brings an error msg
A1spd=spd(wasGoing<0&wasGoing>-pi/4);
A2spd=spd(wasGoing<pi/4&wasGoing>0);
B1spd=spd(wasGoing<pi/2&wasGoing>pi/4);
B2spd=spd(wasGoing<3*pi/4&wasGoing>pi/2);
C1spd=spd(wasGoing<pi&wasGoing>3*pi/4);
C2spd=spd(wasGoing<(-3*pi/4)&wasGoing>-pi);
D1spd=spd(wasGoing<(-pi/2)&wasGoing>(-3*pi/4));
D2spd=spd(wasGoing<(-pi/4)&wasGoing>(-pi/2));
meanA1spd=mean(A1spd);
meanA2spd=mean(A2spd);
meanB1spd=mean(B1spd);
meanB2spd=mean(B2spd);
meanC1spd=mean(C1spd);
meanC2spd=mean(C2spd);
meanD1spd=mean(D1spd);
meanD2spd=mean(D2spd);
y=[meanA2spd,meanB1spd,meanB2spd,meanC1spd,meanC2spd,meanD1spd,meanD2spd,meanA1spd];
figure(10);
bar(y); % or could do hist(speeds) for ugly version
xlabel('run heading (degrees)');
ylabel('run speed (mm/min)');

% NO: displ=eset.gatherSubField('run','pathLength'); %gets all the distances of the runs
% runDur=eset.gatherSubField('run','runTime'); %gets all the run durations in frames
% meanThetas=eset.gatherSubField('run','meanTheta');
% A1displ=displ(prevDir<0 && prevDir>(7*pi)/4);
% A1runDur=runDur(maxAngle<0 && maxAngle>(7*pi)/4);
% A1speed=mean((A1displ*lengthPerPixel)/(A1runDur/240));
% A2displ=displ(maxAngle<pi/4 && maxAngle>0);
% A2runDur=runDur(maxAngle<pi/4 && maxAngle>0);
% A2speed=mean((A2displ*lengthPerPixel)/(A2runDur/240));
% B1displ=displ(maxAngle<pi/2 && maxAngle>pi/4);
% B1runDur=runDur(meanTheta<pi/2 && maxAngle>pi/4);
% B1speed=mean((B1displ*lengthPerPixel)/(B1runDur/240));
% B2displ=displ(maxAngle<(3*pi)/4 && maxAngle>pi/2);
% B2runDur=runDur(maxAngle<(3*pi)/4 && maxAngle>pi/2);
% B2speed=mean((B2displ*lengthPerPixel)/(B2runDur/240));
% C1displ=displ(maxAngle<pi && maxAngle>(3*pi)/4);
% C1runDur=runDur(maxAngle<pi && maxAngle>(3*pi)/4);
% C1speed=mean((C1displ*lengthPerPixel)/(C1runDur/240));
% C2disp=displ(maxAngle<(5*pi)/4 && maxAngle>pi);
% C2runDur=runDur(maxAngle<(5*pi)/4 && maxAngle>pi);
% C2speed=mean((C2displ*lengthPerPixel)/(C2runDur/240));
% D1displ=displ(maxAngle<(3*pi)/2 && maxAngle>(5*pi)/4);
% D1runDur=runDur(maxAngle<(3*pi)/2 && maxAngle>(5*pi)/4);
% D1speed=mean((D1displ*lengthPerPixel)/(D1runDur/240));
% D2displ=displ(maxAngle<(7*pi)/4 && maxAngle>(3*pi)/2);
% D2runDur=runDur(maxAngle<(7*pi)/4 && maxAngle>(3*pi)/2);
% D2speed=mean((D2displ*lengthPerPixel)/(D2runDur/240));
%
%