% Script that analyzes tracks from movies of crawling larvae 
% subjected to spatial temperature gradient.
%
% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) files, is loaded into MATLAB, then processed to show how larval behavior depends on 
% temperature
% 
%% START WITH SOME ANALYSIS PARAMETERS
minTrackLength = 500; % only accept tracks more than X frames in duration.
lasttime = 1200; % time where we stop being interested (in seconds).
cuts = [-inf inf]; % where the offset steps are.
lengthPerPixel = 1; %pixels/sec
%lengthPerPixel = 0.1054; %mm/pixel with the 8X lens (L2/L3); also correct boundary values
lengthPerPixel = 0.0686; %12X lens (L1); also correct boundary values

%% LOAD THE FILE(S)
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles();
end

%% CLEAN THE TRACKS
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

%% FIX HEAD-TAIL ORIENTATION, SET SEGMENT SPEEDS, DO SEGMENTATION (set
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

%% Figure 1 (1B: thermotaxis index for each stage at each temp range)
%Set a 0.4" strip in the center where x=0
positions=eset.gatherField('loc');
xpositions=positions(1,:); 
hotFrames=sum(xpositions<1226); %enter the boundary value for warm frames here
coldFrames=sum(xpositions>1374); %enter the boundary value for cool frames here
thindex=(hotFrames-coldFrames)/(hotFrames+coldFrames) %use in Thindex.m 

%Calculate error bars for the thindex:
numTracks=eset.gatherField('trackNum');
wtIndexTrack=zeros(1,length(numTracks));
prefIndexTrack=zeros(1,length(numTracks));
for i=1:length(numTracks);
    pos=eset.expt.track(i).getDerivedQuantity('loc');
    posX=pos(1,:);
    numXhot=sum(posX<1226); %enter boundary value again
    numXcold=sum(posX>1374); %enter boundary value again
    prefIndexTrack(i)=(numXhot-numXcold)/(numXhot+numXcold); 
    wtIndexTrack(i)=sum(posX<1226|posX>1374); %enter boundary values
end

prefIndexTrack=prefIndexTrack(wtIndexTrack>0);
wtIndexTrack=wtIndexTrack(wtIndexTrack>0);

wtIndexTrackNorm=wtIndexTrack/sum(wtIndexTrack); %normalized
wtMean=sum(wtIndexTrackNorm.*prefIndexTrack); %results in a weighted mean (which matches TI)
wtMeanStd=sqrt(sum(wtIndexTrackNorm.*(prefIndexTrack-wtMean).*(prefIndexTrack-wtMean)));
wtMeanSqrt=sqrt(length(prefIndexTrack)-1);
thindexError=wtMeanStd/wtMeanSqrt 
    
%% Figure 2 (5K: maximum head sweep angle for each direction quadrant)
wasGoing=eset.gatherSubField('headSwing','prevDir');
maxAngle=eset.gatherSubField('headSwing','maxTheta');
A=maxAngle(wasGoing<(pi/4)&wasGoing>(-pi/4));
Amean=mean(abs(A)); 
Ameandeg=radtodeg(Amean); 
B=maxAngle(wasGoing<(3*pi/4)&wasGoing>(pi/4));
Bmean=mean(abs(B));
Bmeandeg=radtodeg(Bmean);
C=maxAngle(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4));
Cmean=mean(abs(C));
Cmeandeg=radtodeg(Cmean);
D=maxAngle(wasGoing<(-pi/4)&wasGoing>(-3*pi/4));
Dmean=mean(abs(D));
Dmeandeg=radtodeg(Dmean);
x=[1,2,3,4];   
y=[Ameandeg,Bmeandeg,Cmeandeg,Dmeandeg];
figure(5);
bar(x,y);
xlabel('direction quadrant (degrees)');
set(gca,'xtickLabel',{'-45 to +45','+45 to +135','+135 to -135','-135 to -45'}); %convention is that 0 degrees is to the right
ylabel('average maximum head sweep angle (degrees)');
title('average maximum head sweep angle by direction quadrant');

%% Figure 3 (5J: probability of starting new run during each head sweep when
% facing towards warmth or cold) 
accs=eset.gatherSubField('headSwing','accepted');
hd=eset.gatherSubField('headSwing','headDir');
Aaccs=accs(hd<(pi/4)&hd>(-pi/4));
Aaccsyes=Aaccs(Aaccs>0); 
Anew=length(Aaccsyes)/length(Aaccs);
Caccs=accs(hd>(3*pi/4)|hd<(-3*pi/4));
Caccsyes=Caccs(Caccs>0);
Cnew=length(Caccsyes)/length(Caccs);
y=[Anew,Cnew];
figure(6);
bar(y);
xlabel('direction quadrant (degrees)');
set(gca,'xtickLabel',{'-45 to +45','+135 to -135'});
ylabel('probability of starting new run during each head sweep (%)');
title('likelihood of accepting head sweep by direction quadrant');

%% Figure 4 (5I: probability that head sweeps first towards warmth) 
firstmax=eset.gatherSubField('firsths','maxTheta');
firstprevDir=eset.gatherSubField('firsths','prevDir');
warmFrames=sum((firstmax>0&firstprevDir<pi)|(firstmax<0&firstprevDir>pi));
coolFrames=sum((firstmax>0&firstprevDir>pi)|(firstmax<0&firstprevDir<pi));
warmFirst=(warmFrames/(warmFrames+coolFrames))*100;
figure(7);
bar(warmFirst);
ylabel('probability head sweeps first towards warm direction (%)');
set(gca,'xtickLabel',{[]});

%% Figure 5 (5H: number of head sweeps during each turning event)
numberHS = eset.gatherSubField('reorientation','numHS');
numberHS1=numberHS(numberHS>0);
figure(8);
hist(numberHS1,[1:10]);
xlabel('number of head sweeps during each turning event');
ylabel('count');

%% Figure 6 (5E: percent of runs in each heading direction (degrees))
aveThetas=eset.gatherSubField('run','meanTheta'); 
A1testing=aveThetas(aveThetas<0&aveThetas>-pi/4);
A2testing=aveThetas(aveThetas<pi/4&aveThetas>0);
B1testing=aveThetas(aveThetas<pi/2&aveThetas>pi/4);
B2testing=aveThetas(aveThetas<3*pi/4&aveThetas>pi/2);
C1testing=aveThetas(aveThetas<pi&aveThetas>3*pi/4);
C2testing=aveThetas(aveThetas<(-3*pi/4)&aveThetas>-pi);
D1testing=aveThetas(aveThetas<(-pi/2)&aveThetas>(-3*pi/4));
D2testing=aveThetas(aveThetas<(-pi/4)&aveThetas>(-pi/2));
numRuns=length(aveThetas);
A1runt=((length(A1testing))/numRuns)*100; 
A2runt=((length(A2testing))/numRuns)*100;
B1runt=((length(B1testing))/numRuns)*100;
B2runt=((length(B2testing))/numRuns)*100;
C1runt=((length(C1testing))/numRuns)*100;
C2runt=((length(C2testing))/numRuns)*100;
D1runt=((length(D1testing))/numRuns)*100;
D2runt=((length(D2testing))/numRuns)*100;
y=[A2runt,B1runt,B2runt,C1runt,C2runt,D1runt,D2runt,A1runt];
figure(9);
bar(y);
xlabel('run heading (degrees)');
set(gca,'xtickLabel',{'-45 to 0','0 to +45','+45 to +90','+90 to +135','+135 to +180','+180 to -135','-135 to -90','-90 to -45'});
ylabel('count (%)');

%% Figure 7 (5D: average run speed(mm/min) for each run heading)
spdraw=eset.gatherField('speed','run'); %units of sdpraw is pixels/s
spd=spdraw*lengthPerPixel*60; %convert to mm/min
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
bar(y);
xlabel('run heading (degrees)');
set(gca,'xtickLabel',{'-45 to 0','0 to +45','+45 to +90','+90 to +135','+135 to +180','+180 to -135','-135 to -90','-90 to -45'});
ylabel('run speed (mm/min)');

%% Figure 8 (5B: run duration (s) for various run headings (degrees))
runDur=eset.gatherSubField('run','runTime');
aveThetasdeg=radtodeg(aveThetas);
figure(11);
plot(aveThetasdeg,runDur,'b*');
xlabel('run heading (degrees)');
ylabel('run duration (s)'); 

%% Figure 9 (5C: run duration (s)) 
ArunDur=runDur(aveThetas<(pi/4)&aveThetas>(-pi/4));
CrunDur=runDur(aveThetas>(3*pi/4)|aveThetas<(-3*pi/4));
figure(12);
subplot(1,2,1), hist(ArunDur); 
xlabel('duration (s)');
ylabel('number of runs in A direction quadrant (toward cool end)');
subplot(1,2,2), hist(CrunDur);
xlabel('duration (s)');
ylabel('number of runs in C direction quadrant (toward warm end)'); 
for i=1:20;
    ArunDurhist(i)=sum((ArunDur>5*i-5)&(ArunDur<5*i));
    ArunDurhistfrac(i)=ArunDurhist(i)/length(ArunDur);
    CrunDurhist(i)=sum((CrunDur>5*i-5)&(CrunDur<5*i));
    CrunDurhistfrac(i)=CrunDurhist(i)/length(CrunDur);
end
figure(15);
xAxis=[5:5:100];
plot(xAxis,ArunDurhistfrac);
xlabel('duration(s)');
ylabel('fraction of runs in A direction quadrant (toward cool end)');
figure(16);
plot(xAxis,CrunDurhistfrac);
xlabel('duration(s)');
ylabel('fraction of runs in C direction quadrant (toward warm end)');

%% Figure 10 (5G: new heading AFTER each turn for each direction
% quadrant)
chg=hd-wasGoing;
for i=1:length(chg);
    if chg(i)>pi;
        chg(i)=chg(i)-2*pi;
    end
    
    if chg(i)<-pi;
        chg(i)=chg(i)+2*pi;
    end
end
Achg=chg(wasGoing<(pi/4)&wasGoing>(-pi/4));
Aaccepteds=accs(wasGoing<(pi/4)&wasGoing>(-pi/4));
AchgAccs=Achg(Aaccepteds>0);
Achgdeg=radtodeg(AchgAccs);
Bchg=chg(wasGoing<(3*pi/4)&wasGoing>(pi/4));
Baccepteds=accs(wasGoing<(3*pi/4)&wasGoing>(pi/4));
BchgAccs=Bchg(Baccepteds>0);
Bchgdeg=radtodeg(BchgAccs);
Cchg=chg(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4)); 
Caccepteds=accs(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4));
CchgAccs=Cchg(Caccepteds>0);
Cchgdeg=radtodeg(CchgAccs); 
Dchg=chg(wasGoing<(-pi/4)&wasGoing>(-3*pi/4));
Daccepteds=accs(wasGoing<(-pi/4)&wasGoing>(-3*pi/4));
DchgAccs=Dchg(Daccepteds>0);
Dchgdeg=radtodeg(DchgAccs);

Awas=wasGoing(wasGoing<(pi/4)&wasGoing>(-pi/4));
AwasAccs=Awas(Aaccepteds>0);
AwasAccsdeg=radtodeg(AwasAccs);
Bwas=wasGoing(wasGoing<(3*pi/4)&wasGoing>(pi/4));
BwasAccs=Bwas(Baccepteds>0);
BwasAccsdeg=radtodeg(BwasAccs);
Cwas=wasGoing(wasGoing>(3*pi/4)|wasGoing<(-3*pi/4));
CwasAccs=Cwas(Caccepteds>0);
CwasAccsdeg=radtodeg(CwasAccs);
Dwas=wasGoing(wasGoing<(-pi/4)&wasGoing>(-3*pi/4));
DwasAccs=Dwas(Daccepteds>0);
DwasAccsdeg=radtodeg(DwasAccs);
figure(13);
subplot(1,4,1), plot(Achgdeg,AwasAccsdeg,'b.'); 
xlabel('heading change after each turn (degrees)');
ylabel('heading before each turn (degrees) (A)');
subplot(1,4,2), plot(Bchgdeg,BwasAccsdeg,'b.');
xlabel('heading change after each turn (degrees)');
ylabel('heading before each turn (degrees) (B)');
subplot(1,4,3), plot(Cchgdeg,CwasAccsdeg,'b.');
xlabel('heading change after each turn (degrees)');
ylabel('heading before each turn (degrees) (C)');
subplot(1,4,4), plot(Dchgdeg,DwasAccsdeg,'b.');
xlabel('heading change after each turn (degrees)');
ylabel('heading before each turn (degrees) (D)');

%% Figure 11 (5F: heading change by the end of each run for each direction
% quadrant)
startAng=eset.gatherSubField('run','startTheta');
endAng=eset.gatherSubField('run','endTheta');
delt=(endAng-startAng);
for i=1:length(delt);
    if delt(i)>pi;
        delt(i)=delt(i)-2*pi;
    end
    
    if delt(i)<-pi;
        delt(i)=delt(i)+2*pi;
    end
end
Adelt=delt(aveThetas<(pi/4)&aveThetas>(-pi/4));
Adeltdeg=radtodeg(Adelt);
Bdelt=delt(aveThetas<(3*pi/4)&aveThetas>(pi/4));
Bdeltdeg=radtodeg(Bdelt);
Cdelt=delt(aveThetas>(3*pi/4)|aveThetas<(-3*pi/4));
Cdeltdeg=radtodeg(Cdelt);
Ddelt=delt(aveThetas<(-pi/4)&aveThetas>(-3*pi/4));
Ddeltdeg=radtodeg(Ddelt);
Asub=aveThetas(aveThetas<(pi/4)&aveThetas>(-pi/4));
Asubdeg=radtodeg(Asub);
Bsub=aveThetas(aveThetas<(3*pi/4)&aveThetas>(pi/4));
Bsubdeg=radtodeg(Bsub);
Csub=aveThetas(aveThetas>(3*pi/4)|aveThetas<(-3*pi/4));
Csubdeg=radtodeg(Csub);
Dsub=aveThetas(aveThetas<(-pi/4)&aveThetas>(-3*pi/4));
Dsubdeg=radtodeg(Dsub);
figure(14);
subplot(1,4,1), plot(Adeltdeg,Asubdeg,'b.');
xlabel('heading change by the end of each run (degrees)');
ylabel('heading before each turn (degrees) (A)');
subplot(1,4,2), plot(Bdeltdeg,Bsubdeg,'b.');
xlabel('heading change by the end of each run (degrees)');
ylabel('heading before each turn (degrees) (B)');
subplot(1,4,3), plot(Cdeltdeg,Csubdeg,'b.');
xlabel('heading change by the end of each run (degrees)');
ylabel('heading before each turn (degrees) (C)');
subplot(1,4,4), plot(Ddeltdeg,Dsubdeg,'b.');
xlabel('heading change by the end of each run (degrees)');
ylabel('heading before each turn (degrees) (D)');
