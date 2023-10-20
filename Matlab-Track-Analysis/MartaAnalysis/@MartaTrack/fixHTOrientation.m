function fixHTOrientation(mt, varargin)
%makes sure each head and tail point is aligned
%by aligning only the tails (less likely to move)
%function fixHTOrientation(track, varargin)

detectInvalid = true;
debug = false;
varargin = assignApplicable(varargin);

pt = mt.pt;
inds = find([pt.htValid]);
oldspine = pt(inds(1)).spine;

for j = 2:length(inds)
    
    sp = pt(inds(j)).spine;
    d1 = sum(sqrt(sum((sp - oldspine).^2)));
    d2 = sum(sqrt(sum((sp(:,end:-1:1) - oldspine).^2)));
    if (d2 < d1)
        sp = sp(:,end:-1:1);
        pt(inds(j)).spine = sp;
     
        pt(inds(j)).head = sp(:,1);
        pt(inds(j)).tail = sp(:,end);
    end
    oldspine = sp;
 end

vel = diff([pt(inds).mid], [], 2);
mh = [pt(inds(2:end)).head] - [pt(inds(2:end)).mid];
midtail = [pt(inds(2:end)).tail] - [pt(inds(2:end)).mid];
if (sum(dot(vel, mh)) < sum(dot(vel,midtail)))
    for j = 1:length(pt)
        temp = pt(j).head;
        pt(j).head = pt(j).tail;
        pt(j).tail = temp;
        sp = pt(j).spine;
        sp = sp(:,end:-1:1);
        pt(j).spine = sp;
    end
end
%pt(1)
mt.pt = pt;
mt.dq = [];

%now look for points where the spine jumps a signficant amount compared to
%the contour length, & mark as invalid
if (detectInvalid)
    mt.markHTInvalid(0.05);
    mt.fixHTOrientation('detectInvalid', false, varargin{:});
    return;
end

%now let's segment track into regions of continuously valid ht & consider
%those individually
if (any (~[pt.htValid]))
    htv = [pt.htValid];
    start = [1 find(diff(htv) > 0)];
    stop = [find(diff(htv) < 0) - 1, length(htv)];
    for j = 1:length(start)
        inds = start(j):stop(j);
        if (length(inds) < 6) %why bother
            continue;
        end
        vel = diff([pt(inds).mid], [], 2);
        mh = [pt(inds(2:end)).head] - [pt(inds(2:end)).mid];
        midtail = [pt(inds(2:end)).tail] - [pt(inds(2:end)).mid];
        
        if (sum(dot(vel, mh)) < sum(dot(vel,midtail)))
            for k = inds
                
                temp = pt(k).head;
                pt(k).head = pt(k).tail;
                pt(k).tail = temp;
                sp = pt(k).spine;
                sp = sp(:,end:-1:1);
                pt(k).spine = sp;
            end
        end
    end
    mt.pt = pt;
end

%on last run through, use median filtering to fix parts that continue to be
%bad


mintime = 5; %seconds minimum time for something to be bad before it needs to be fixed
track = mt;


if (sum([pt.htValid]) < 4 || sum([pt.htValid]) < 4*mintime/track.dr.interpTime)
    return; %track is too short
end
sp = track.getDerivedQuantity('smoothSpeed');
speedthresh = percentile(sp, 0.33);
dpthresh = -0.5;
assignApplicable(varargin);
minpts = mintime/track.dr.interpTime;
dp = track.getDerivedQuantity('vel_dp');

inds = find(sp > speedthresh & isfinite(dp));

if (length(inds) < 2)
    disp ('no points above speed cut in fixHTOrientation');
    return;
end

dp = double(interp1(inds,dp(inds), 0:(length(dp)+1), 'nearest', 'extrap'));
dpfilt = medfilt1(dp, minpts);
badinds = dpfilt < dpthresh;
badinds = find(imdilate(badinds, ones([1 3])));
badinds = badinds((badinds >= 1) & (badinds <= length(track.dq.eti))); 
badtimes = track.dq.eti((badinds)); %time when it's going backwards

%find the original points that are in that time range
%note we are assuming that the interpolation time <= sampling time, which
%should be valid; we space out by thirds to make sure we cover all points
%note unique sorts, so there is no need to do that here
badtimes = unique([badtimes, badtimes-track.dr.interpTime/3, badtimes+track.dr.interpTime/3]);
badinds = unique(interp1([pt.et], 1:track.npts, badtimes,'nearest', 'extrap'));

%swap head and tail in bad time range
if (~isempty(badinds))
    for j = badinds
        temp = pt(j).head;
        pt(j).head = pt(j).tail;
        pt(j).tail = temp;
        sp = pt(j).spine;
        sp = sp(:,end:-1:1);
        pt(j).spine = sp;      
    end
    mt.pt = pt;
end