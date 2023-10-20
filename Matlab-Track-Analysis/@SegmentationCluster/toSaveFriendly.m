function sc2 = toSaveFriendly(sc)

if (length(sc) > 1)
    for j = 1:length(sc);
        sc2(j) = toSaveFriendly(sc(j));
    end
    return;
end

sc2 = SegmentationCluster();
f = fieldnames(sc);
f = setdiff(f, {'knownPoints'});



for j = 1:length(f)
    sc2.(f{j}) = sc.(f{j});
end

kp = sc.knownPoints;
if (isempty(kp) || isempty(kp(1).track) || isempty(kp(1).inds))
    np.filename = '';
    np.locInFile = [];
    sc2.knownPoints = np;
    return;
end
t = [kp.track];
e = unique([t.expt]);

for j = 1:length(e)
    np(j).filename = e(j).fname;
    tr = t([t.expt] == e(j));
    np(j).locInFile = uint32([]);
    for k = 1:length(tr);
        inds = [kp([kp.track] == tr(k)).inds];
        pt = tr(k).pt(unique(tr(k).getDerivedQuantity('mapinterpedtopts', false, inds)));
        np(j).locInFile = [np(j).locInFile uint32([pt.locInFile])];
    end
end

sc2.knownPoints = np;