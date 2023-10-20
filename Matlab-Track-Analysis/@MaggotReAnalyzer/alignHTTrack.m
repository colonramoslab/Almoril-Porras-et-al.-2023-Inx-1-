function alignHTTrack (mra, track)
%function alignHTTrack (mra, track)

existsAndDefault('track', mra.track);
valid = double([track.pt.htValid]);

start = find(diff(valid) > 0) + 1;
if (valid(1))
    start = [1 start];
end

stop = find(diff(valid) < 0);
if (valid(end))
    stop = [stop length(valid)];
end
pt = track.pt;
for j = 2:length(pt)
    if (mra.alignHTPt(pt(j), pt(j-1)))
        temp = pt(j).head;
        pt(j).head = pt(j).tail;
        pt(j).tail = temp;
    end
end
track.calculateDerivedQuantity({'iloc', 'sloc'}, true);
track.calculateDerivedQuantity('vel', true);

vi = track.getDerivedQuantity('mapPtsToInterped');
if (mra.debug)
    figure(10); clf(10);
    track.calculateDerivedQuantity({'ihead', 'itail', 'imid', 'ibodytheta'}, true);
    mh = track.dq.ihead - track.dq.imid;
    et = track.getDerivedQuantity('eti');
    th = track.getDerivedQuantity('theta');
    bth = atan2(mh(2,:), mh(1,:));
    plot (et, th, 'b-', et(find(valid)), th(find(valid)), 'b.', et(find(valid)), bth(find(valid)), 'r.');
    title ('heading and body angle vs. time pre-align');
end

for j = 1:length(start)
    inds = start(j):stop(j);
    if (length(inds) < 4)
        continue;
    end
    v = track.dq.vel(:,vi(inds));
    %size(inds)
    th = [pt(inds).head] - [pt(inds).tail];
    norm_th = sqrt(sum(th.^2, 1));
    th = th ./ [norm_th;norm_th];

    dp = dot(v,th);

    if (sum(dp) < 0)
        for k = inds
            temp = pt(k).head;
            pt(k).head = pt(k).tail;
            pt(k).tail = temp;
        end
    end
end
track.pt = pt;
if (mra.debug)
    figure(11); clf(11);
    track.calculateDerivedQuantity({'ihead', 'itail', 'imid', 'ibodytheta'}, true);
    mh = track.dq.ihead - track.dq.imid;
    et = track.getDerivedQuantity('eti');
    th = track.getDerivedQuantity('theta');
    bth = atan2(mh(2,:), mh(1,:));
    plot (et, th, 'b-', et(logical(valid)), th(logical(valid)), 'b.', et(logical(valid)), bth(logical(valid)), 'r.');
    title ('heading and body angle vs. time post-align');
end