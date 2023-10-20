function condenseKnownPoints(sc)
%packs all known points from a single track into a single known points
%structure
%function condenseKnownPoints(sc);

if (length(sc) > 1)
    for j = 1:length(sc)
        sc(j).condenseKnownPoints();
    end
    return
end

kp = sc.knownPoints;
if (isempty(kp))
    return;
end
tr = unique([kp.track]);

for j = 1:length(tr)
    kpnew(j).track = tr(j);
    kpnew(j).inds = sort([kp([kp.track] == tr(j)).inds]);
end
    
sc.knownPoints = kpnew;

end

