function setSegmentSpeeds (expt, mso)
% function setSegmentSpeeds (track, mso)

% If other segmentation options (mso) are given, use those.  Otherwise just
% default to existing track.so options.  
existsAndDefault('mso', expt.so);
expt.so = mso;

% Extract all the speeds and curvatures from the track (or is it all
% tracks?)
sp = expt.gatherField(mso.speed_field);
cv = expt.gatherField('curv');

% Keep only the curvatures that are above the segmentation threshold.
highcurv = abs(cv) > mso.curv_cut;
% Issue warnings if there are zero or few points above threshold.  
if (isempty(highcurv))
    disp(['locInFile = ' num2str(track.locInFile) ' no high curvature']);
    return
end
if (sum(highcurv) < 4)
    disp(['locInFile = ' num2str(track.locInFile) ' few high curvature points']);
end

% Calculate average speed (and standard deviation) of larvae with
% curvatures above the threshold.  Then set the "stop" speed cutoff to be
% the mean + stdev.  
u = mean(sp(highcurv));
s = std(sp(highcurv));
expt.so.stop_speed_cut = u + s;

% Not sure what this does.  Expands which points are above the curvature
% threshold, then cuts out the original points above threshold.  The result
% being a set of points that occur near when the larval curvature is above
% threshold?  Then take the average speed (and stdev) of these points,
% using that as the cutoff for restarting a run?  
nearhc = imdilate(highcurv, ones(5)) &~highcurv;
u = mean(sp(nearhc));
s = std(sp(nearhc));
expt.so.start_speed_cut = u + s;

[expt.track.so] = deal(expt.so);


