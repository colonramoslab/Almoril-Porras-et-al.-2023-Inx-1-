%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   TEMPORAL LIGHT GRADIENT ANALYSIS SCRIPT
%
%   Written by Liz    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to analyze larval tracks navigating a temporal light gradient
% produced by a digital projector.

% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and light (.bin) files, is loaded into
% MATLAB, then processed to show how larval behavior depends on 
% temporal light gradients.

% Set minimum and maximum light threshold values using target light levels
maxLT = 200;    %this should be projector target light level max that was used to aquire data
minLT = 0;      %this should be projector target light level min that was used to aquire data 
onThresh = minLT + .95*(maxLT-minLT);   %max thresh is 5% less of diff between max/min
offThresh = minLT + .05*(maxLT-minLT);  %min thresh is 5% more of diff between max/min
endSec = 1400; %cut data 23 minutes (derived by looking at all time and deciding that response was max during first 23 mins)
trimBuffer=50; %pixel buffer distance to trim within extraction window

%create resetAll as a var to reload a new experiment from scratch

%% LOAD EXPERIMENT
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
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP AND SEGMENT TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %the below code can be used to display the eset cleaning levels 
% %graphically and to change default values specified below
% ecl.getReport(eset);

% %throw out first minute of data acquistion
% eset.executeExperimentFunction('trimTracks',[60 1200],[]);

%Only look at first 1400 secs of data - determined this was the range
%where larvae remain maximally sensitive to light changes
%have to do this before you do anything 
%trim a little past 1400 so that reorientation plots for last time bin
%aren't erroneous

%analysis options struct keeps track of which analyses have been run on the
%data you have open

if ~exist('ao','var')
    ao=[];
end

if ~bypass
   
    if ~isfield(ao,'trim') || ao.trim==false
        disp(['trimming track time to be less than ' num2str(endSec) ' and smaller than extraction window']);
        loc=eset.gatherSubField('pt','loc'); %gather all center locations
        minPair=min(loc,[],2)';
        maxPair=max(loc,[],2)';
        eset.executeExperimentFunction('trimTracks',[0 endSec+100],[minPair+trimBuffer maxPair-trimBuffer]);
        ao.trim = true;
    end

    if ~isfield(ao,'clean') || ao.clean==false
        ecl = ESetCleaner;
        ecl.minHTValid= .98;
        ecl.minDist=75;
        ecl.minSpeed = 0.75;
        ecl.minPts = 500;
        
        disp('cleaning eset');
        ecl.clean(eset);
        ao.clean=true;
    end
    
    %close figures opened by ecl.clean(eset)
    close all

    %write to file
    writeL=input('Do you want to write to mat file? 1 or 0: ');
    if writeL
        if~exist('toFile','var')|| exist('reSave','var') && reSave ||  exist('resetAll','var')
            disp('writing to mat file');
            fname=eset.expt(1).fname;
            inds=strfind(fname, '\');
            fstub=[fname(1:inds(end)) 'MatFiles\'];
            mkdir(fstub);
            eset.toMatFiles([fstub 'Temporal_trim_' num2str(endSec) '_clean_' num2str(trimBuffer)]);
        end
    end
end

% %fixes head/tail orientation
% 
% if ~isfield(ao,'fixht') || ao.fixht==false
%     disp('fixing head/tail orientation');
%     eset.executeTrackFunction('fixHTOrientation');
%     ao.fixht=true;
% end

%set headsweep theta min to 15 deg for segmentation (found this doesn't
%miss as many headsweeps for second instars
for j=1:length(eset)
    eset.expt(j).so.headswing_start=deg2rad(15);
end

% Set segmentation threshold speeds
if ~isfield(ao,'setSS') || ao.setSS==false
    disp('setting segment speeds');
    eset.executeTrackFunction('setSegmentSpeeds');
    ao.setSS=true;
end

% Segment the tracks
if ~isfield(ao,'seg') || ao.seg==false
    disp('segmenting tracks');
    eset.executeTrackFunction('segmentTrack');
    ao.seg=true;
end

%Calc Gradient Type
%Add light on/off information
eti=[eset(1).expt(1).globalQuantity(1).xData];    %gather time info
lt=[eset(1).expt(1).globalQuantity(1).yData]; 
dlt=[eset(1).expt(1).globalQuantity(2).yData];    %gather derivative target light info

%calculate conversion factor from frames to minutes
secPerFrame=eset(1).expt(1).dr.interpTime;
frame2min=(60/secPerFrame);

% Determine gradient type
derZero=find(dlt==0); %find indicies where derivative of target light was zero (there will be many of these for square waves)
if length(derZero) > 20
    gradientType = 'Square';
else
    gradientType = 'Triangle';
end

%Add light on/off counting vectors (useful for looking at beh metric for
%only a certain time period after light value has changed)
switch (gradientType)
    case ('Square')
        if ~isfield(eset.expt(1).track(1).dq, 'periodCount') || isfield(ao,'resetAll')
            disp ('deriving light time information');
            binLT = lt > 1; %thresh light values
            onInds = find(diff(binLT)>0); %find inds where light goes on
            offInds= find(diff(binLT)<0); %find inds where light goes off
            % plot(eti,binLT,eti(onInds),binLT(onInds),'r.',eti(offInds),binLT(offInds),'g.')
            
            %make vectors of how long has passed since light was on/off
            %by convention, the period starts when the lights go on, so
            %periodCount is the same as onTimeCount
            periodCount=ones(size(binLT));
            offTimeCount=ones(size(binLT));
            
            periodCount(onInds)=0;
            offTimeCount(offInds)=0;
            
            %fill in lights on/off vectors
            count=1;
            for j=1:length(periodCount)-1
                if periodCount(j)==1
                    periodCount(j)=count;
                    count=count+secPerFrame;
                elseif periodCount(j)==0
                    count=1;
                end
            end
            
            count=1;
            for j=1:length(offTimeCount)-1
                if offTimeCount(j)==1
                    offTimeCount(j)=count;
                    count=count+secPerFrame;
                elseif offTimeCount(j)==0
                    count=1;
                end
            end
            
            %determine period
            period=round(max(periodCount));
            
            %make first half period invalid by substiting in large
            %negative numbers
            periodCount(1:period/2)=-100;
            offTimeCount(1:period/2)=-100;
            
            %Get rid of any half periods at the end
            endpt=max(onInds)+1;
            periodCount(endpt:length(periodCount))=-100;
            offTimeCount(endpt:length(offTimeCount))=-100;
            
            %Add these vectors as a global quantity
            
            disp ('adding light global quantity');
            for j = 1:length(eset.expt) %add quantity to all expts
                gq = GlobalQuantity();
                
                %add lights on time count
                gq.xField = 'eti'; % label
                gq.xData = eti;     % elapsed time, defined above
                
                gq.fieldname = 'periodCount'; % label
                gq.yData = periodCount;        % time elapsed since the lights have gone on in sec
                eset.expt(j).addGlobalQuantity(gq); %add global quantity
                
                gq.fieldname = 'tLightOff'; % label
                gq.yData = offTimeCount;         % time elapsed since the lights have gone off in sec
                eset.expt(j).addGlobalQuantity(gq);
            end
        end
    case('Triangle')
        if ~isfield(eset.expt(1).track(1).dq, 'periodCount') || isfield(ao, 'resetAll')
            disp ('deriving light time information');
            bindlt=dlt>0;
            onInds = find(diff(bindlt)>0); %find inds where light starts increasing
            offInds= find(diff(bindlt)<0); %find inds where light starts decreasing
            
            %Make period counter vector
            periodCount=ones(size(lt));
            periodCount(onInds)=0;
            
            %fill in period counter vector
            count=1;
            for j=1:length(periodCount)-1
                if periodCount(j)==1
                    periodCount(j)=count;
                    count=count+secPerFrame;
                elseif periodCount(j)==0
                    count=1;
                end
            end
            
            %determine period
            period=ceil(max(periodCount));
            
            %Get rid of any half periods at the end by making large
            %negative
            % endpt=max(onInds)+1;
            %periodCount(endpt:length(periodCount))=-100;
            
            %get rid of first period since animals were just put down
            %periodCount(1:period)=-100;
            
            %Add these vectors as a global quantity
            disp ('adding light global quantity');
            for j = 1:length(eset.expt) %add quantity to all expts
                gq = GlobalQuantity();
                
                %add period count count
                gq.xField = 'eti'; % label
                gq.xData = eti;     % elapsed time, defined above
                
                gq.fieldname = 'periodCount'; % label
                gq.yData = periodCount;        % time elapsed since the lights have gone on in sec
                eset.expt(j).addGlobalQuantity(gq); %add global quantity
            end
        end
end

% Need to check segmentation by playing some movies here
if ~isfield(ao,'ao.checkSeg');
    checkSeg=input('Do you want to check segmentation? Enter 1 or 0: ');
    if checkSeg
        %use below for watching movies in order of ihtV score
        %figure;
        %eset.makeHistogram('ihtValid',0.8:0.01:1,'mean');
        %ihtV = eset.gatherField('ihtValid','mean');
        %[~,I] = sort(ihtV);
        %figure;
        for j=1:length(eset.expt)
            for i=1:length(eset.expt(j).track)
                eset.expt(j).track(i).playMovie;
                pause
            end
        end
    end
    ao.checkSeg=true;
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize fignum to zero
if checkSeg
    close all %close segmentation checking figures
end

fignum = 0;

%Calculate beh metrics for intervals around light changes
dOnThresh=40;
dOffThresh=-dOnThresh;
secInt=20;
secIntToFrame=20*(1/secPerFrame);

%Determine figure axes
binSize=period/10;
ptx=0+binSize/2:binSize:period; %xperiod axis
etx=0+binSize/2:binSize:endSec; %time axis x

% % Plot light value and derivative of light value
% fignum = fignum + 1; figure(fignum); clf(fignum);
% plot(eti,lt,eti,dlt, '-r');
% axis([0 endSec -20 205])
% xlabel('eti');
% legend('Target Light', 'dTarget Light');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT TURNING/PAUSING RATE AS A FUNCTION OF TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot TURNING RATE V. TIME WITH LIGHT
fignum = fignum + 1; figure(fignum); clf(fignum);
[trh,trh_eb]=eset.makeReorientationHistogram('eti',etx,'minHS',1,'incllastrun',true); %include last run since will be same for temporal assay
ax = plotyy(etx, trh, eti, lt);
hold on
%errorbar(etx,trh,trh_eb);
xlim (ax(1), [0 endSec]); %set both x axes to be proper range
xlim (ax(2), [0 endSec]);
xlabel('time')
axes(ax(1)); ylabel('turning rate 1/minutes'); %set left axis name
axes(ax(2)); ylabel('target light value'); %set right axis name
title('turning');
hold off

%Plot turning rate V. Period Count fignum=fignum+1; figure(fignum); clf(fignum);
[a,eb]=eset.makeReorientationHistogram('periodCount',ptx,'minHS',1,'incllastrun',true);
%gather light and time vectors that correspond to the tLightOn
%range by gathering time info for one period from the first on time
%to the second on time with a five point buffer
xp=eti(onInds(1)-5:onInds(2)+5);
%shift times back to zero so axis can start at zero
xp=xp-xp(1);
%get light values that correspond to the x time
l=lt(onInds(1)-5:onInds(2)+5);
[ax,h1,h2]=plotyy(ptx,a,xp,l);
hold on;
errorbar(ptx,a,eb);
set(ax(1),'xlim',[0 period]); %set both x axes to be proper range
set(ax(2),'xlim',[0 period]);
ylabel(ax(1),'turning rate 1/minutes'); %set left axis name
ylabel(ax(2),'target light value'); %set right axis name
xlabel('time (period count)');
title('turning');
hold off

%Plot PAUSING RATE V. TIME WITH LIGHT
fignum = fignum + 1; figure(fignum); clf(fignum);
hold on
[trh,eb]=eset.makeReorientationHistogram('eti',etx,'maxHS',0,'incllastrun',true); %include last run since will be same for temporal assay
ax = plotyy(etx, trh, eti, lt); 
%errorbar(etx,trh,eb);
xlim (ax(1), [0 endSec]); %set both x axes to be proper range
xlim (ax(2), [0 endSec]);
xlabel('time')
ylabel(ax(1),'pausing rate 1/minutes'); %set left axis name
ylabel(ax(2),'target light value'); %set right axis name
title('pausing');
hold off

%Plot PAUSING RATE V. TIME SINCE LIGHT HAS BEEN ON
fignum=fignum+1; figure(fignum); clf(fignum);
[a,eb]=eset.makeReorientationHistogram('periodCount',ptx,'maxHS',0, 'inclastrun',true);
ax=plotyy(ptx,a,xp,l);
hold on
errorbar(ptx,a,eb);
xlim (ax(1), [0 period]); %set both x axes to be proper range
xlim (ax(2), [0 period]);
xlabel('time (period count)');
ylabel(ax(1),'pausing rate 1/minutes'); %set left axis name
ylabel(ax(2),'target light value'); %set right axis name
title('pausing');
hold off


% %plot PAUSING RATE AND TURNING RATE TOGETHER V. TIME SINCE LIGHT HAS BEEN
% %ON
% h1 = eset.makeReorientationHistogram('tLightOn', 0:5:period, 'minHS',1);
% h2 = eset.makeReorientationHistogram('tLightOn', 0:5:period, 'maxHS',0);
% fignum=fignum+1; figure(fignum); clf(fignum);
% plot (0:5:period, h1, 0:5:period, h2);

%Look at speed v. time since light on
fignum=fignum+1; figure(fignum); clf(fignum);
eset.meanField2vsField1('periodCount','speed',0:binSize:period);
xlabel('time (period count)');
ylabel('speed');
title('speed v. time');


%look at number of head sweeps v. time since light on
reoPer=eset.gatherFromSubField('reorientation','periodCount','position','start');
numHS=eset.gatherSubField('reorientation','numHS');
[x,meany,stderr,stddev]=meanyvsx(reoPer(numHS>0), numHS(numHS>0), 0:binSize:period);
fignum=fignum+1; figure(fignum); clf(fignum);
%plotyy(x,meany,xp,l);
%xlim (ax(1), [0 period]); %set both x axes to be proper range
%xlim (ax(2), [0 period]);
%xlabel('Time since light turned on');
%axes(ax(1)); ylabel('Mean number of head sweps per reorientation'); %set left axis name
%axes(ax(2)); ylabel('target light value'); %set right axis name
plot(x,meany);
hold on
errorbar(x,meany,stderr);
xlim([0 period]);
ylabel('Mean number of headsweeps per reorientation');
xlabel('time (period count)');
hold off

%look at head sweep death v. time since light on
hsStartTimes=eset.gatherFromSubField('headSwing','periodCount','position','start');
hsAng=eset.gatherSubField('headSwing','maxTheta');
[x,meany,stderr,stddev]=meanyvsx(hsStartTimes, rad2deg(abs(hsAng)), 0:binSize:period);
fignum=fignum+1; figure(fignum); clf(fignum);
plot(x,meany);
hold on
errorbar(x,meany,stderr);
xlabel('time (period count)');
ylabel('Mean max headsweep angle');
hold off



 