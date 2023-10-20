% Script that analyzes tracks from movies of crawling larvae 
% subjected to temporally varying (and spatially uniform) temperature.
% Used to be called "analyzeMason_".  
%
% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and temperature (.tmp) files, is loaded into
% MATLAB, then processed to show how larval behavior depends on 
% temperature.
% 
% 
% ANALYSIS PARAMETERS 
%  
minTrackLength = 500; % only accept tracks more than X frames in duration.
lasttime = 2400; % time where we stop being interested (in seconds).
cuts = [-inf inf]; % where the offset steps are.
tempSensorNum = 7; % this starts at zero, if you look at a .tmp file.
lengthPerPixel = 1;
%lengthPerPixel = 0.118; % (with the 8X lens)
lengthPerPixel = 0.0766; % (with the 12X lens)
%lengthPerPixel = ? % (with the even zoomier lens, not used yet)
%
% dT/dt must be >/< than +/- threshold for events to be counted as occurring
% during rising/falling temperatures.  
dtempthresh = 1E-3; %degree/sec


% LOAD THE FILE(S)
%
if (~exist('eset','var'))
    eset=ExperimentSet.fromFiles();
end

% CLEAN THE TRACKS
%
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    ecl = ESetCleaner;
    ecl.minHTValid = 0.90;
    ecl.minDist = 50;
    ecl.minSpeed = 0.5;
    ecl.minPts = minTrackLength;
    ecl.clean(eset);
    
    cleanEset = false;
end

% LOAD AND ADD TEMPERATURE INFORMATION

% Load temperature/times from the .tmp file(s), store in cells
if (~exist('temps','var'))
      
   for j = 1:length(eset.expt)
       %
       % Find the .tmp file name (fn2), same name as the .bin file.  
       fn = eset.expt.fname;
       ind = strfind(fn,'.');
       fn2 = [fn(1:ind) 'tmp'];
       %
       % Load the values contained in the .tmp file
       tempdata = load(fn2);
       %
       % Extract times and temperatures from the right rows in the matrix.   
       times{j} = tempdata(tempSensorNum*3+2,:)/1000; % convert to seconds
       temps{j} = tempdata(tempSensorNum*3+3,:);
   end
   clear tempdata;

end

% Convert tracks to the TemperatureMaggotTrack class.    
for j = 1:length(eset.expt)
    
    clear t2;     
    if (~isa(eset.expt(j).track(1), 'TemperatureMaggotTrack'))
        for k = 1:length(eset.expt(j).track)
            t2(k) = TemperatureMaggotTrack(eset.expt(j).track(k));
        end
        eset.expt(j).track = t2; 
    end
       
end

% Add temperature to the eset as a global quantity 
if (~isfield(eset.expt(1).track(1).dq, 'temp'))
   disp ('assigning temperature information');
   
   for j = 1:length(eset.expt)
       
       % 't' is a vector holding all the unique values of elapsed time.
       % 'tm' is a vector holding all the corresponding temperatures.
       % (there are many redundant points in the .tmp files) 
       [t,I] = unique(times{j} - times{j}(1));
       tm = temps{j}(I);
              
       %make an interpolated time that is of the finest resolution of the
       %time sampling
       ti = min(t):min(abs(diff(t))):max(t);
       
       %too much information?  no interpolated temperature
       tmi = interp1(t, tm, ti);
     
       % smooth the data with a Gaussian filter
       lowpasstime = 10;
       deltatime = diff(ti(1:2)); 
       sigma = lowpasstime / deltatime;
       tms = lowpass1D(tmi, sigma);
       dtm = deriv(tms, 1); %since the data is already smoothed we can 
                           % always use a value of 1 here
       
       % Add temperature (T) and time derivative of temperature (dT/dt)
       % as global quantities.
       gq = GlobalQuantity();
       
       gq.xField = 'eti'; % label
       gq.xData = ti;     % elapsed time, defined above
       %default: gq.derivationMethod = @GlobalQuantity.oneDinterpolation;
       
       gq.fieldname = 'temp'; % label
       gq.yData = tms;        % temperature, defined above
       eset.expt(j).addGlobalQuantity(gq);
       
       % Replace 'temp' with 'dtemp' and add that as a gq also.  
       gq.fieldname = 'dtemp'; % label
       gq.yData = dtm;         % derivative of temperature, defined above
       eset.expt(j).addGlobalQuantity(gq);
       % Note that this addGlobalQuantity() function adds 'temp' and
       % 'dtemp' as expt.track(_).dq.temp in each track, as well as being
       % global quantities of expt( ).  

   end
   
end


% FIX HEAD-TAIL ORIENTATION, CORRECT SPEED, DO SEGMENTATION

% Head-tail orientation fixing.  
existsAndDefault('fixht','true');
if (fixht)
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

% Speed correction (both temperature and larval size) 
% note: currently a bug where if you run this twice it messes up the size
% correction somehow.  
if (~exist('aspeed', 'var') || (exist('respeed','var') && respeed)) 
    
    figure(1); clf(1);
    figure(2); clf(2);
    
    for j = 1:length(eset.expt)
        
        % A: SPEED VS. TEMPERATURE CORRECTION
        figure(1); clf(1);
        
        % Grab speed from every point of every track, and the temperature
        % of the plate corresponding to these speeds.  Store as vectors.
        temp = eset.expt(j).gatherField('temp');
        speed = eset.expt(j).gatherField('speed');
        
        % Linear fit to a speed vs. temperature plot
        p = polyfit(temp, speed, 1);
        % tempToSA (SA=speed adjust?) is a vector that covers the range of
        % temperatures in the experiment, in 1 deg. C increments.  
        clear tempToSA;
        tempToSA(1,:) = floor(min(temp)-0.01):1:ceil(max(temp)+0.01);
        % Speeds from the speed vs. temperature linear fit that correspond
        % to the temperatures in tempToSA.  Then a new row in tempToSA is
        % the average speed divided by the sfit values.
        sfit = polyval(p, tempToSA(1,:));
        tempToSA(2,:) = mean(speed)./sfit;

        % Plot the full speed vs. temperature values, then also the linear
        % fit (in integer steps of temperature).  
        plot (temp, speed, 'k.', 'MarkerSize', 1); hold on
        plot (tempToSA(1,:), sfit,'r-','LineWidth', 3);
        xlabel('temp');
        ylabel('speed px/sec');
        
        
        % Store tempToSA as a property of the track
        [eset.expt(j).track.tempToSA] = deal(tempToSA);
        
        % Using calculateDerivedQuantity (for @TemperatureMaggotTrack) to 
        % find adjusted speed.  Uses calculateAdjustedSpeed, which uses
        % adj. speed = speed x (avg. speed)/(speed at temp. in sfit).  
        eset.expt(j).executeTrackFunction('calculateDerivedQuantity', 'adjusted_speed');
        
        % B. SPEED VS. LARVA SIZE CORRECTION
        
        figure(2); clf(2);
        
        % Obtain the area (in pixels?) of the larvae.
        % [evaluateTrackExpression goes through each track and executes 
        % whatever string is in its argument, and the resulting value (or 
        % vector) is stored]. 
        % In this particular case, looks like each track is given a vector
        % of size [size(track.dq.eti)] (i.e., number of frames) filled with
        % the average area of the larva over the course of the whole track.
        area = eset.expt(j).evaluateTrackExpression('ones(size(track.dq.eti))*mean([track.pt.area])');
        % Obtain the adjusted speed from every point in every track, the 
        % same as was done above for regular speed.  Note that area and
        % aspeed are vectors with the same dimensions, it's just that area
        % has been assigned the avg. area for the whole track, not allowed
        % to vary from point to point.  
        aspeed = eset.expt(j).gatherField('adjusted_speed');
        % Extract average area for each track, and the average adjusted 
        % speed for each track.  I don't think these are used below, so we
        % can probably eliminate this.  (?)
        meanarea = eset.expt(j).evaluateTrackExpression('mean([track.pt.area])');
        meanspeed = eset.expt(j).evaluateTrackExpression('mean(track.getDerivedQuantity(''adjusted_speed''))');

        % Linear fit to adjusted speed vs. larval area.  
        p2 = polyfit(area, aspeed, 1);
        sizeToSA(1,:) = floor(min(area)-0.01):1:ceil(max(area)+0.01);
        szfit = polyval(p2, sizeToSA(1,:));
      
        % and put in a figure...
        plot (area, aspeed, 'k.', 'MarkerSize', 1); hold on
        plot (sizeToSA(1,:), szfit,'r-','LineWidth', 3);
        xlabel('area of maggot px^2');
        ylabel('speed px/sec');
        
        % The average speed (already adjusted for temperature) of all
        % larvae in all tracks.  
        mas = mean(aspeed);
        
        % cf (correction factor?) is the average speed divided by the speed
        % from the (speed vs. area) fit, at the area in the specific track.
        % 
        % Then adjust the adjusted_speed by multiplying by the cf of each
        % track.  
        for k = 1:length(eset.expt(j).track)
            cf = mas/polyval(p, mean([eset.expt(j).track(k).pt.area]));
            eset.expt(j).track(k).dq.adjusted_speed = eset.expt(j).track(k).dq.adjusted_speed*cf;
        end
        
        % Use the adjusted speed as the speed field in segmentation options
        % (so, see @MaggotSegmentOptions).  
        eset.expt(j).so.speed_field = 'adjusted_speed';
        % Change the curvature threshold (for deciding when a run ends)
        % Not sure if need this!
        eset.expt(j).so.curv_cut = 1;
        %
        % Rewrites the experiment segmentation options (so) onto the track
        % segmentation options.  What for?
        [eset.expt(j).track.so] = deal(eset.expt(j).so);
        
    end
    
    respeed = false;
    
end

% Set segmentation threshold speeds (uses the new adjusted speed)
existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% Segment the tracks
existsAndDefault('segment', true);
if (segment)
    eset.executeTrackFunction('segmentTrack');
    segment = false; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DONE WITH EXTRACTION.  NOW MAKE ALL THE FIGURES WE WANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLOT TEMPERATURE AND DERIVATIVE OF TEMPERATURE
figure(3); clf(3);

subplot(2,1,1);
% Temperature.  Include the raw and smoothed versions.  
hold on
plot(ti,tmi);
plot(ti,tms,'r-','LineWidth',2);
xlabel('time (sec)');
ylabel('agar temperature (C)');
title('Temperature!');
legend('raw temperature','smoothed temperature');
hold off

% dT/dt.  Include +/- the heating/cooling threshold.  
subplot(2,1,2);
dtempthreshVec = ones(1,length(ti))*dtempthresh;
plot(ti,dtm,ti,dtempthreshVec,ti,-dtempthreshVec);
xlabel('time (sec)');
ylabel('d/dt of temperature (C/sec)');
title('derivative of temperature!');



% PLOT REORIENTATION RATES DURING HEATING AND COOLING
%
reo=eset.gatherField('reorientation');
reoHS = [reo.numHS];
%reo = [reo(reoHS>0)];
clear reo;
%
reodT = eset.gatherFromSubField('reorientation','dtemp','position','mean');
reoT = eset.gatherFromSubField('reorientation','temp','position','mean');
reodT = reodT(reoHS>0);
reoT = reoT(reoHS>0);
%
alldT = eset.gatherField('dtemp','run','notlast');
allT = eset.gatherField('temp','run','notlast');
%
%
figure(4); clf(4);
%
% Step through the temperature offset sections
for j = 1:(length(cuts) - 1)
    %
    % count number of reorientations during heating/cooling, divide by 
    % number of frames during heating/cooling within a specified
    % temperature range.  Multiply by 240 to get turns/minute.  
    %
    heatRate(j) = 240 * ( sum(reodT>dtempthresh & reoT>cuts(j) & reoT<cuts(j+1))/sum(alldT>dtempthresh & allT>cuts(j) & allT<cuts(j+1)) );
    coolRate(j) = 240 * ( sum(reodT<-dtempthresh & reoT>cuts(j) & reoT<cuts(j+1))/sum(alldT<-dtempthresh & allT>cuts(j) & allT<cuts(j+1)) );
    
    % average temperature between the specified range.  
    tempOffsets(j) = mean(allT(allT > cuts(j) & allT < cuts(j+1)));
    
    % average heating/cooling rates between the specified range (deg C/min)
    dtempHeatRate(j) = 60 * mean(alldT(alldT>dtempthresh & allT>cuts(j) & allT<cuts(j+1)));
    dtempCoolRate(j) = 60 * mean(alldT(alldT<-dtempthresh & allT>cuts(j) & allT<cuts(j+1)));
    
end

% Plot the turning rates during rising/falling temperatures, for each
% offset temperature.  
plot (tempOffsets, heatRate, 'r.-', tempOffsets, coolRate, 'b.-','LineWidth',2,'MarkerSize',30);
xlabel ('<temperature>');
ylabel ('reorientation rate (min^{-1})');
title ('reorientation rate vs. temperature');
legend('heating', 'cooling');
%embiggen();  

% (NEW) alternate method to extracting turning rates.  This finds the
% turning rate for each track, so that we end up with a distribution.  
% numTracks = length(eset.expt.track);
% trackturnsHeat = zeros(1,numTracks);
% trackturnsCool = zeros(1,numTracks);
% %trackTempPref = zeros(1,numTracks);
% trackframesHeat = zeros(1,numTracks);
% trackframesCool = zeros(1,numTracks);
% trackLengths = zeros(1,numTracks);
% weightTempPref = zeros(1,numTracks);
% 
% 
% % First, find the number of heating/cooling turns in each track
% 
% for j=1 : length(reo)
%     
%     if(reodT(j)>dtempthresh)
%         num = reo(j).track.trackNum;
%         trackturnsHeat(num) = trackturnsHeat(num) + 1;
%     end
%     
%     if(reodT(j)<-dtempthresh)
%         num = reo(j).track.trackNum;
%         trackturnsCool(num) = trackturnsCool(num) + 1;
%     end
%     
% end
% 
% % Second, find the number of track points during heating/cooling (and
% % during runs) in each track.
% 
% for j=1 : numTracks
%    
%     for k=1 : eset.expt.track(j).npts
%         
%         if(eset.expt.track(j).isrun(k)>0)
%             
%             if(eset.expt.track(j).dq.dtemp(k)>dtempthresh)
%                 trackframesHeat(j) = trackFramesHeat(j)+1;
%             end
%             
%             if(eset.expt.track(j).dq.dtemp(k)<-dtempthresh)
%                 trackframesCool(j) = trackFramesCool(j)+1;
%             end 
%                      
%         end
%                 
%     end
%     
% end
% 
% % Third, do some more stuff
% 
% 
% 
% % for j=1 : numTracks
% %     trackTempPref(j) = ( trackturnsHeat(j) - trackturnsCool(j) )/( trackturnsHeat(j) + trackturnsCool(j));
% %     trackLengths(j) = eset.expt.track(j).npts;
% %     weightTempPref(j) = trackTempPref(j)*trackLengths(j);
% % end
% % 
% % weightTempPref = weightTempPref / sum(trackLengths);
% % meanTempPref = mean(trackTempPref)
% % meanTempPref2 = mean(weightTempPref)
% % stdevTempPref = std(trackTempPref)
% % sterrTempPref = stdevTempPref / sqrt(numTracks)



% PLOT HEAD SWEEP DEPTH FOR HEATING AND COOLING

figure(5); clf(5);

% angle of the first head sweep of each reorientation
maxtheta = eset.gatherSubField('firsths', 'maxTheta');
thetaBins = 0:30:180;


for j=1:(length(cuts)-1)
    
    % sweep depth during heating and cooling (absolute value, in degrees)
    maxthetaHeat = rad2deg(abs(maxtheta(reodT>dtempthresh & reoT>cuts(j) & reoT<cuts(j+1))));
    maxthetaCool = rad2deg(abs(maxtheta(reodT<-dtempthresh & reoT>cuts(j) & reoT<cuts(j+1))));
    
    maxthetaHeatAVG(j) = mean(maxthetaHeat);
    maxthetaCoolAVG(j) = mean(maxthetaCool);
        
end

plot(tempOffsets, maxthetaHeatAVG, 'r.-',tempOffsets, maxthetaCoolAVG, 'b.-','LineWidth',2,'MarkerSize',30);
xlabel('<temperature>');
ylabel('AVG of first HS angle (deg)');
title('Head sweep angle');
legend('heating','cooling');


% histograms of sweep depth, for only the first temperature offset
maxthetaHeat = rad2deg(abs(maxtheta(reodT>dtempthresh & reoT>cuts(1) & reoT<cuts(2))));
maxthetaCool = rad2deg(abs(maxtheta(reodT<-dtempthresh & reoT>cuts(1) & reoT<cuts(2))));
maxthetaHeatHist = hist(maxthetaHeat,thetaBins);
maxthetaCoolHist = hist(maxthetaCool,thetaBins);

% Plot these heating and cooling histograms
figure(6); clf(6);
plot(thetaBins, maxthetaHeatHist/sum(maxthetaHeatHist), 'r-', thetaBins, maxthetaCoolHist/sum(maxthetaCoolHist), 'b-','Linewidth', 2);
heatLegend = ['heating (avg = ' num2str(maxthetaHeatAVG(1),2) ' deg)'];
coolLegend = ['cooling (avg = ' num2str(maxthetaCoolAVG(1),2) ' deg)'];
legend(heatLegend,coolLegend);
xlabel('head sweep angle (deg)');
ylabel('number of head sweeps (norm.)');
title('Head sweep distribution at first temp. offset');


% FIND NUMBER OF HEAD SWINGS FOR HEATING AND COOLING

figure(7); clf(7);
reoHS2 = reoHS(reoHS>0);

for j=1:(length(cuts)-1)
    
    reoHSheat = reoHS2(reodT>dtempthresh & reoT>cuts(j) & reoT<cuts(j+1));
    reoHScool = reoHS2(reodT<-dtempthresh & reoT>cuts(j) & reoT<cuts(j+1));
    % 
    % Average number of head sweeps per reorientation during heating/cooling.
    reoHSheatMean(j) = mean(reoHSheat);
    reoHScoolMean(j) = mean(reoHScool);

    % Percentage of reorientations with #HS>1, during heating/cooling
    reoHSheatPercent(j) = 100*sum(reoHSheat>1)/length(reoHSheat);
    reoHScoolPercent(j) = 100*sum(reoHScool>1)/length(reoHScool);
         
end

plot(tempOffsets, reoHSheatMean, 'r.-',tempOffsets, reoHScoolMean, 'b.-','LineWidth',2,'MarkerSize',30);
xlabel('<temperature>');
ylabel('AVG #HS per reorientation');
title('Number of head sweeps');
legend('heating','cooling');

figure(8); clf(8);
plot(tempOffsets, reoHSheatPercent, 'r.-',tempOffsets, reoHScoolPercent, 'b.-','LineWidth',2,'MarkerSize',30);
xlabel('<temperature>');
ylabel('% of reorientations with #HS>1');
title('Reorientations with # head sweeps > 1');
legend('heating','cooling');


% PRINT (LATER WRITE TO FILE) ALL INFORMATION WE WANT TO KEEP

% number of tracks (also want # of animals, maybe not from Matlab):
numTracks = length(eset.expt.track);
display(['number of tracks: ' num2str(numTracks)]);

% average speed
meanSpeed = mean(speed);
meanSpeedUnits = meanSpeed * lengthPerPixel;
display(['average speed: ' num2str(meanSpeed,3) ' px/s (' num2str(meanSpeedUnits,3) ' mm/s)']);
for j=1:(length(cuts)-1)
    speedOffset = speed(temp>cuts(j) & temp<cuts(j+1));
    meanSpeedOffset(j) = mean(speedOffset);
    meanSpeedOffsetUnits(j) = lengthPerPixel * meanSpeedOffset(j);
end
display(['average speed (at each offset): ' num2str(meanSpeedOffset,3) ' px/s (' num2str(meanSpeedOffsetUnits,3) ' mm/s)']);


% average temperature(s), and heating/cooling rate(s)
display(['temperature offset(s): ' num2str(tempOffsets,3) ' C']);
display(['  heating rate(s):  ' num2str(dtempHeatRate,3) ' C/min']);
display(['  cooling rate(s): ' num2str(dtempCoolRate,3) ' C/min']);

% reorientation rate(s) (heating/cooling)
display('reorientation rate per larva');
display(['  heating: ' num2str(heatRate,3) ' turns/min']);
display(['  cooling: ' num2str(coolRate,3) ' turns/min']);

% average first HS sweep depth (heating/cooling)
display('average first head sweep depth');
display(['  heating: ' num2str(maxthetaHeatAVG,3) ' deg']);
display(['  cooling: ' num2str(maxthetaCoolAVG,3) ' deg']);

% average # HS per reorientation (heating/cooling)
display('average number of head sweeps per reorientation');
display(['  heating: ' num2str(reoHSheatMean,3)]);
display(['  cooling: ' num2str(reoHScoolMean,3)]);

% % of reorientations with more than 1 HS (heating/cooling)
display('Percent of reorientations with more than one head sweep');
display(['  heating: ' num2str(reoHSheatPercent,3)]);
display(['  cooling: ' num2str(reoHScoolPercent,3)]);

% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIGURE!
% figure(9); clf(9)
% 
% t = [eset.expt.track];
% hs = [t.headSwing];
% hsaccepted = [hs.accepted];
% hstemp = eset.gatherField('temp', 'headSwing', 'mean');
% hsdtemp = eset.gatherField('dtemp', 'headSwing', 'mean');
% 
% for j = 1:(length(cuts) - 1)
%     hsrising = (hstemp > cuts(j) & hstemp < cuts(j+1) & hsdtemp > dtempthresh);
%     hsfalling = (hstemp > cuts(j) & hstemp < cuts(j+1) & hsdtemp < dtempthresh);
%     risingaccepted(j) = mean(hsaccepted(hsrising));
%     fallingaccepted(j) = mean(hsaccepted(hsfalling));
%     hsmeantemp(j) = mean(hstemp(hsrising | hsfalling));
% end
% 
% bar (hsmeantemp', [fallingaccepted;risingaccepted]');
% set(gca, 'XTick', hsmeantemp);
% title ('probability of accepting a headsweep');
% legend('falling','rising')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%warning(warning_state);

return



%figure(3)
%bar ((uptemp+downtemp)/2, downfrac./upfrac-1);
%set(gca, 'XTick', (uptemp+downtemp)/2);
%xlabel ('<temperature>');
%ylabel ('headsweep ratio - 1');
%title ('Fraction of time headsweeping when dtemp < 0 / Fraction of time headsweeping when dtemp > 0 - 1');
%embiggen();