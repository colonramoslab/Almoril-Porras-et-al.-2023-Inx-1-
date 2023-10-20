% version 2 of analyzeMason.m.  
%
% Script that analyzes tracks from movies of crawling larvae 
% subjected to temporally varying (and spatially uniform) temperature.  
%
% Data from track extraction software (.bin files), along with their 
% associated timing (.tim) and temperature (.tmp) files, is loaded into
% MATLAB, then processed to show how larval behavior depends on 
% temperature
% 
basedir = 'J:\BEHAVIOR\05252010\Analysis\';
% Directory where all the .bin, .tim, .tmp files are located.  Note the 
% need for the semicolon at the end of these lines, and for the last 
% backslash.  Doesn't work without this.  Later we probably want to 
% replace this part with something that just pops up and asks the user 
% for a location.
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Load Data from Files %%%%%%
%
% d is an array of file names in basedir, specifically files of .bin format.
% Entries in d include .name, .date, .bytes, .isdir, .datenum.  
% Here we use .name to define the path to the file, and .bytes to display
% the file sizes
d = dir([basedir '*.bin']);
disp (['total size on disk is ' num2str(sum([d.bytes]))]);
% Choose which file in the basedir folder to load.  This can be a single
% number ('1' would load the first .bin file and associated .tim and .tmp
% files -- 'first' meaning first alphabetically, or date, or however they
% are sorted), or a vector specifying which files (each of which becomes
% an experiment) we want to load.  
whichToLoad = 1;
% Don't load tracks that are less than 500 frames?
% Or does this mean don't load whole experiments less than 500 frames?
minptstoload = 500;
% Keep warnings from filling the screen if something goes wrong?
warning_state = warning('off', 'all');

% LOAD TRACKS AND TIMING
% Step through each experiment we want to load... then load it!
% (Skip if any experiments were already loaded -- i.e., the 'expt'
%  variable already exists)
if (~exist('expt','var'))
   for j = 1:length(whichToLoad)
       memory % display memory available
       tic
       % Obtain the .bin file name (fn) and the .tim file name (fn2).
       % Assumed that these files have the same names.
       fn = [basedir d(whichToLoad(j)).name];
       ind = strfind(fn,'.');
       fn2 = [fn(1:ind(end)) 'tim'];
       % Command that loads the experimental data from the .bin file.
       % fromFile is a static method of the Experiment class.  
       expt(j) = Experiment.fromFile(fn, fn2, false, [], minptstoload); 
       toc % display time it took to load track data
   end
end

% LOAD TEMPERATURE INFORMATION
% Step through each experiment and load temperature data from the .tmp
% file.
% (Skip if this was already done -- i.e., 'temps' exists already)
if (~exist('temps','var'))
      
   for j = 1:length(whichToLoad)
       
       % Find the .tmp file name (fn2), assumed to have the same name
       % as the .bin file.  
       fn = [basedir d(whichToLoad(j)).name];
       ind = strfind(fn,'.');
       fn2 = [fn(1:ind) 'tmp'];
       % Load the values contained in the .tmp file
       % (might want to rename this -- 'data' is kind of generic
       data = load(fn2);
       % Extract times and temperatures from the appropriate
       % columns in the matrix.  These numbers might vary depending on 
       % which temperature sensor we want values from.  
       % Also, should we clear 'data' from memory afterward?
       times{j} = data(end-4,:)/1000; % convert to s from ms?
       temps{j} = data(end-3,:);
       
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Eliminate slow and/or short tracks %%%%%%

% Minimum average speed of a track allowed, and the minimum rms 
% displacement allowed.  Tracks not meeting the criteria are eliminated.
speedthresh = 1; %pixels/sec
distthresh = 75; %pixels excursion

% Only run if 'retrim' does not exist yet, or if it has been reset to true.  
if (~exist('retrim', 'var') || retrim)
    
   for j = 1:length(expt)
       
       % Calculate average speed.  'meanspeed' is a vector, which contains
       % the average speed for each track.  
       meanspeed = expt(j).evaluateTrackExpression('mean(track.getDerivedQuantity(''speed''))');
       % Only keep tracks whose average speed is above threshold
       expt(j).track = expt(j).track(meanspeed > speedthresh); %#ok<SAGROW>
       
       % Calculate the rms displacement.  'dt' (distance travled?) is a
       % vector holding the rms displacement for each track.
       dt = sqrt(expt(j).evaluateTrackExpression('max(sum(track.getDerivedQuantity(''displacement'').^2))'));
       % Only keep tracks whose rms displacment is above threshold.
       expt(j).track = expt(j).track(dt > distthresh); %#ok<SAGROW>
       
       clear t2; % wouldn't we want to also clear this afterward?
       % If the first track (and by implication, all tracks in the expt)
       % does not belong to the TemperatureMaggotTrack class, then
       % convert it to this class.  Done by setting each element of 't2'
       % to the TemperatureMaggotTrack version of the track, then 
       % replacing the original tracks with 't2'.  
       if (~isa(expt(j).track(1), 'TemperatureMaggotTrack'))
           for k = 1:length(expt(j).track)
               t2(k) = TemperatureMaggotTrack(expt(j).track(k));
           end
           expt(j).track = t2; %#ok<SAGROW>
       end
       
   end
   retrim = false;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Misc.  %%%%%%

% If there is no experiment set defined yet, then make one.
if (~exist('eset','var'))
    
   eset = ExperimentSet();
   eset.expt = expt; 
   % should we clear 'expt' after this?  Is expt taking up double memory?
   
   % Run the head-tail orientation routine.  I think this goes through
   % all the tracks.  Whenever the dot product of the orientation and
   % the velocity is negative (or more negative than some threshold) AND
   % the speed of the animal is above some other threshold, the orientation
   % is wrong, and will be changed.  
   disp('fixing HT orientation - this can take a while');
   eset.executeTrackFunction('fixHTOrientation');
   %
   %eset.executeExperimentFunction('stitchTracks', 5, 15);
   % Why is this commented out now?  Seems like a standard step.  
   
end

% Incorporate the temperatures extracted above from the .tmp files into
% the eset.expt, as a global quantity.  
% 
% Run if 'temp' is not already a field in .dq (derived quantity).  I think
% this part might need to change now, since temperatures are not put there
% anymore.  
if (~isfield(expt(1).track(1).dq, 'temp'))
   disp ('assigning temperature information; this should be reasonably fast');
   
   for j = 1:length(eset.expt)
       
       % Define 't' as a vector holding all the unique values of elapsed
       % time (not sure why there would be duplicates).  Then define 'tm'
       % as a vector holding all the temperatures corresponding to those
       % elapsed times.  
       
       %find all the unique times the temperature was taken
       [t,I] = unique(times{j} - times{j}(1));
       
       %grab the temperature at those times
       tm = temps{j}(I);
       
       %make an interpolated time that is of the finest resolution of the
       %time sampling
       ti = min(t):min(abs(diff(t))):max(t);
       
       %too much information?  no interpolated temperature
       tmi = interp(t, tm, ti);
       
       
       % Both functions lowpass1D and deriv are in the 'basic routines'
       % folder.  I assume this is mostly about smoothing data.  
       
       %now let's smooth the data with a gaussian filter of width sigma,
       %which we'll set in seconds
       lowpasstime = 10;
       deltatime = diff(ti(1:2)); 
       sigma = lowpasstime / deltatime;
       tms = lowpass1D(tmi, sigma);
       dtm = deriv(tm, 1); %since the data is already smoothed we can always use a value of 1 here
       
       
       % Adds temperature and time derivative of temperature as global
       % quantities in the experiment set.  A global quantity is a new
       % field that is known based on the time (and maybe position) of a
       % track.  xField is a label, and xData holds the values.  In this
       % case, xField/Data is elapsed time, and fieldname/yData is
       % temperature.  
       gq = GlobalQuantity();
       
       gq.xField = 'eti'; % label
       gq.xData = ti;      % elapsed time, defined above
       %default: gq.derivationMethod = @GlobalQuantity.oneDinterpolation;
       
       gq.fieldname = 'temp'; % label
       gq.yData = tms;         % temperature, defined above
       eset.expt(j).addGlobalQuantity(gq);
       
       % Replace 'temp' with 'dtemp' and add that as a gq also.  
       gq.fieldname = 'dtemp'; % label
       gq.yData = dtm;         % derivative of temperature, defined above
       eset.expt(j).addGlobalQuantity(gq);
       % Note that this addGlobalQuantity() function basically performs
       % gq.addQuantityToTrack(expt.track(j)) for all the tracks.  I guess
       % this makes temperature and the associated time part of the
       % information contained within a track.  In this case, 'temp' and
       % 'dtemp' are added as expt.track(_).dq.temp.  

   end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART IV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Correct for Speed and Larva Size  %%%%%%
%
% Only run if adjusted speed variable doesn't exist, or another variable
% respeed has been set to true.  
if (~exist('aspeed', 'var') || (exist('respeed','var') && respeed)) 
    
    % FIGURE!
    figure(1); clf(1);
    
    for j = 1:length(eset.expt)
        
        % A: SPEED VS. TEMPERATURE CORRECTION
        
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
        
        % Actual correction is here.  
        % Not sure what this first part is for.  Storing tempToSA as a
        % property of the track?  
        [eset.expt(j).track.tempToSA] = deal(tempToSA);
        % Using calculateDerivedQuantity (here the one for
        % @TemperatureMaggotTrack, which is what we have here) to find the
        % adjusted speed.  Uses the function calculateAdjustedSpeed, which
        % is adj. speed = speed x (avg. speed)/(speed at temp. in sfit).  
        eset.expt(j).executeTrackFunction('calculateDerivedQuantity', 'adjusted_speed');
        
        % B. SPEED VS. LARVA SIZE CORRECTION
        
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
        p = polyfit(area, aspeed, 1);
        % The average speed (already adjusted for temperature) of all
        % larvae in all tracks.  
        mas = mean(aspeed);
        
        % Additional adjustment to speed based on the size of the larvae,
        % sort of a shorthand version of what was done above for
        % temperature.  
        %
        % cf (correction factor?) is the average speed divided by the speed
        % from the (speed vs. area) fit, at the area in the specific track.
        % 
        % Then adjust the adjusted_speed by multiplying by the cf of each
        % track.  
        for k = 1:length(eset.expt(j).track)
            cf = mas/polyval(p, mean([expt(j).track(k).pt.area]));
            eset.expt(j).track(k).dq.adjusted_speed = eset.expt(j).track(k).dq.adjusted_speed*cf;
        end
        % Now the speed has been adjusted to account for the fact that
        % larvae move faster at higher temperature and for the fact that
        % bigger larvae move faster.   
        
    end
    
    respeed = false;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Segmentation  %%%%%%

for j = 1:length(eset.expt)
    
    % Do segmentation if track.run does not exist yet, or if the variable
    % resegment exists and is true.
    if (isempty(eset.expt(j).track(1).run) || (exist('resegment','var') && resegment))
        
        % Use the adjusted speed as the speed field in segmentation options
        % (so, see @MaggotSegmentOptions).  
        eset.expt(j).so.speed_field = 'adjusted_speed';
        % Change the curvature threshold (for deciding when a run ends)
        eset.expt(j).so.curv_cut = 1;
        %
        % Rewrites the experiment segmentation options (so) onto the track
        % segmentation options.  What for?
        [eset.expt(j).track.so] = deal(eset.expt(j).so);
        
        % Sets speed thresholds: how slow a larva should move to constitute
        % the end of a run, and how fast they must move to constitute the
        % start of a run.  
        eset.expt(j).executeTrackFunction('setSegmentSpeeds')
        
        % The actual track segmentation, breaking tracks down into runs and
        % reorientations.  These are both functions for the class
        % @MaggotTrack.  
        eset.expt(j).executeTrackFunction('segmentTrack');
        
    end
  
end
resegment = false;

%%%%%%%%%%%%%%%%% MORE PARTS %%%%%%%%%%%%%%%%%%%%%%

% This copies the experiment(s), along with all the associated fields, then
% only keeps valid tracks.  Things that make tracks invalid are the larva
% in a track being too slow, too little time spent in runs (too much in
% reorienatations), and the average run time being too short (compared to 
% all runs in the experiment).  
expt2 = eset.expt;
for j = 1:length(expt2)
    expt2(j) = Experiment();
    fn = fieldnames(expt2);
    for k = 1:length(fn)
        expt2(j).(fn{k}) = eset.expt(j).(fn{k});
    end
    valid = true(size(expt2(j).track));
 %   valid(expt2(j).detectPossibleSegmentationProblems) = false;
 %   expt2(j).track = expt2.track(valid);
end
eset.expt = expt2;
% Should clear expt2 afterward?  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The elapsed time (in seconds) where we stop paying attention to the movie
% frames.  This should of course go near the top of the script, as it will
% depend on the individual movie.  
lasttime = 2400;

% Extracts T, dT/dt, t that occur at the ends of runs (i.e., just prior to
% a turning event).  Only include events occurring before lasttime.  
runendtemp = eset.gatherField('temp','runend');
runenddtemp = eset.gatherField('dtemp', 'runend');
runendtime = eset.gatherField('eti', 'runend');
runendtemp = runendtemp(runendtime < lasttime);
runenddtemp = runenddtemp(runendtime < lasttime);

% Extracts T, dT/dt, t that occur *during* all runs (from all experiments
% in the set, all tracks in each experiment).  Only include points
% occurring before lasttime, as above.  
allruntemp = eset.gatherField('temp','run');
allrundtemp = eset.gatherField('dtemp', 'run');
eti = eset.gatherField('eti','run');
allruntemp = allruntemp(eti < lasttime);
allrundtemp = allrundtemp(eti < lasttime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE!
figure(2); clf(2);

% dT/dt must be >/< |threshold| for events to be counted as occurring
% during rising/falling temperatures.  
dtempthresh = 1E-3; %degree/sec

% Cutoffs (in temperature) between stages of offsets.  Needs to start and
% end with -Inf and Inf I guess.  The other numbers indicate breaks between
% stages of oscillation.  This will be problematic if the offset difference
% is small and the lower end of oscillations at a higher temperature
% overlaps with the high ene of oscillations at a lower temperature.  This
% should, of course, be listed at the top of the script, as it will change
% almost every time.  
cuts =  [-Inf 18.5 20.6 22.4 Inf];

% Step through the temperature offset sections
for j = 1:(length(cuts) - 1)
    
    % Count the number of turning events during rising, and the total 
    % number of points during rising.  
    endrising = sum(runendtemp > cuts(j) & runendtemp < cuts(j+1) & runenddtemp > dtempthresh);
    allrising = sum(allruntemp > cuts(j) & allruntemp < cuts(j+1) & allrundtemp > dtempthresh);
    % Find the average temperature (rt) during the offset section (should
    % put this at the beginning of the for loop).  
    rt(j) = mean(allruntemp(allruntemp > cuts(j) & allruntemp < cuts(j+1)));
    % Turning rate while temperature is rising.  endrising/allrising gives
    % the turning rate in turns per frame.  240 converts this to turns per
    % minute (60 sec/min x 4 frames/sec).  
    risingrate(j) = endrising/allrising * 240;
    
    % Same steps as above, but for falling temperatures.  
    endfalling = sum(runendtemp > cuts(j) & runendtemp < cuts(j+1) & runenddtemp < -dtempthresh);
    allfalling = sum(allruntemp > cuts(j) & allruntemp < cuts(j+1) & allrundtemp < -dtempthresh);
    fallingrate(j) = endfalling/allfalling * 240;
    
end

% Plot the turning rates during rising/falling temperatures, for each
% offset temperature.  
plot (rt, risingrate, 'r.-', rt, fallingrate, 'b.-','LineWidth',3,'MarkerSize',50);
xlabel ('<temperature>');
ylabel ('reorientation rate (min^{-1})');
title ('reorientation rate vs. temperature');
legend('temperature increasing', 'temperature decreasing');
embiggen(); % this just makes the font bigger?  (size 16) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE!
figure(3); clf(3)
% Extract the body angle of the larvae in the experiment set.  
th = eset.gatherField('sbodytheta');
% Subset of the body angles, only including angles above pi/4
turning = abs(th) > pi/4;
% Extract T, dT/dt, and t.  Note that this eti will overwrite the one
% obtained in the figure above.  Might want to give this a different name. 
temp = eset.gatherField('temp');
dtemp = eset.gatherField('dtemp');
eti = eset.gatherField('eti');

% Uses the same "cuts" as above.  
% upinds and downinds are indices for frames where the temperature is
% rising/falling (and when time is prior to lasttime; those 2400s should be
% given that variable name here).  The function "find" returns the indices
% of nonzero elements.  
for j = 1:(length(cuts) - 1)
   upinds{j} = find(eti < 2400 & temp > cuts(j) & temp < cuts(j+1) & dtemp > dtempthresh);
   downinds{j} = find(eti < 2400 & temp > cuts(j) & temp < cuts(j+1) & dtemp < -dtempthresh);
end

% Step through the offset steps.  
for j = 1:length(upinds)
    
   % Average temperature at this offset while the temperature is
   % rising/falling.  These should be almost the same, but who knows.  
   uptemp(j) = mean(temp(upinds{j}));
   downtemp(j) = mean(temp(downinds{j}));
   
   % Number of frames spent with larval body angle > pi/4 divided by the
   % total number of frames (while either rising or falling).  
   upfrac(j) = sum(turning(upinds{j}))/length(upinds{j});
   downfrac(j) = sum(turning(downinds{j})) / length(downinds{j});
   
end

% Plot the fraction of the time spent with body angle > pi/4, comparing
% situations of rising and falling temperatures.  Horiztonal axis will just
% be the offset temperature steps.  
plot (uptemp, upfrac, 'r.-', downtemp, downfrac, 'b.-','LineWidth',3,'MarkerSize',50);
xlabel ('<temperature>');
ylabel ('fraction of time spent headsweeping');
title ('fraction of time when body bend angle > 45 vs. temperature');
legend('temperature increasing', 'temperature decreasing');
embiggen(); % font size increase to 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIGURE!
figure(4); %clf(4)?

t = [eset.expt.track];
hs = [t.headSwing];
hsaccepted = [hs.accepted];
hstemp = eset.gatherField('temp', 'headSwing', 'mean');
hsdtemp = eset.gatherField('dtemp', 'headSwing', 'mean');

for j = 1:(length(cuts) - 1)
    hsrising = (hstemp > cuts(j) & hstemp < cuts(j+1) & hsdtemp > dtempthresh);
    hsfalling = (hstemp > cuts(j) & hstemp < cuts(j+1) & hsdtemp < dtempthresh);
    risingaccepted(j) = mean(hsaccepted(hsrising));
    fallingaccepted(j) = mean(hsaccepted(hsfalling));
    hsmeantemp(j) = mean(hstemp(hsrising | hsfalling));
end

bar (hsmeantemp', [fallingaccepted;risingaccepted]');
set(gca, 'XTick', hsmeantemp);
title ('probability of accepting a headsweep');
legend('falling','rising')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning(warning_state);

return



%figure(3)
%bar ((uptemp+downtemp)/2, downfrac./upfrac-1);
%set(gca, 'XTick', (uptemp+downtemp)/2);
%xlabel ('<temperature>');
%ylabel ('headsweep ratio - 1');
%title ('Fraction of time headsweeping when dtemp < 0 / Fraction of time headsweeping when dtemp > 0 - 1');
%embiggen();