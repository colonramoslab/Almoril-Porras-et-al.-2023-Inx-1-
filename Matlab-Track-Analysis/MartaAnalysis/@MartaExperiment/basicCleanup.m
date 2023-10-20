function basicCleanup(expt, mintime)
% trims out bad sections of track & fixes HT direction to match velocity
% function basicCleanup(expt, mintime)

existsAndDefault('mintime', 5);

if (length(expt) > 1)
    for j = 1:length(expt)
        expt(j).basicCleanup(mintime);
    end
    return;
end


t = [expt.track];
for j = 1:length(t)
    pt = [t(j).pt];
    area = mean([pt.area]);
end
maxarea = mean(area) + 3*std(area);
%maxarea = mean(area(area < maxarea)) + 5*std(area(area < maxarea));

for k = 1:length(expt)
    track = expt(k).track;
    valid = true(size(track));
    et = zeros(size(track));
    for j = 1:length(track)
        valid(j) = track(j).removeCollisionPoints(maxarea);
        et(j) = track(j).pt(end).et - track(j).pt(1).et;
    end
    track = track(valid & et > mintime);
    for j = 1:length(track)
        track(j).fixHTOrientation;
    end
    expt(k).track = track;
 %   expt(k).executeTrackFunction('markHTInvalid', 0.1, 'debug', true);
        
end
    


        

