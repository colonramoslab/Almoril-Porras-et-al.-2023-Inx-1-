function arr = toJavaArr(sc)
%converts segmentation cluster or group of segmentation clusters to 
%a string representation that can be saved to disk
%function [str, hashmap] = toString(sc)

arr = javaArray('java.util.HashMap', length(sc));

for j = 1:length(sc)
    no = convertSc(sc);
    arr(j) = objToHashmap(no);
end

%converts segmentation cluster to text friendly version by
%indentifying known points by file name and location in file, rather
%than by handle and index, which can change
function newobj = convertSc(sc)

f = fieldnames(sc);
f = setdiff(f, {'knownPoints', 'clustCov'});

for j = 1:length(f)
    newobj.(f{j}) = sc.(f{j});
end
newobj.clustCov = sc.clustCov(:);

kp = sc.knownPoints;
t = [kp.track];
e = unique(t.expt);

for j = 1:length(e)
    np(j).filename = e(j).fname;
    tr = t([t.expt] == e(j));
    for k = 1:length(tr);
        inds = [kp([kp.track] == tr(k)).inds];
        pt = tr(k).pt(unique(tr(k).getDerivedQuantity('mapinterpedtopts', inds)));
        np(j).locInFile = [np(j).locInFile [pt.locInFile]];
    end
end

newobj.knownPoints = np;
