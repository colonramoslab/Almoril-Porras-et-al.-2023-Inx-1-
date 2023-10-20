function alignHTSegment (mra, track, inds)
%function alignHTSegment (mra, track, inds)
%
%aligns each point to the previous then aligns the entire segment

track.calculateDerivedQuantity('vel', false);
vi = track.getDerivedQuantity('mapPtsToInterped');
v = track.dq.vel(:,vi(inds));

for j = inds(2:end)
    if (mra.alignHTPt(track.pt(j), track.pt(j-1)))
        temp = track.pt(j).head;
        track.pt(j).head = track.pt(j).tail;
        track.pt(j).tail = temp;
    end
end

th = [track.pt(inds).head] - [track.pt(inds).tail];
norm_th = sqrt(sum(th.^2, 1));
th = th ./ [norm_th;norm_th];
size(v)
size(th)
dp = dot(v,th);

if (sum(dp) < 0)
    for j = inds
        temp = track.pt(j).head;
        track.pt(j).head = track.pt(j).tail;
        track.pt(j).tail = temp;
    end
end

