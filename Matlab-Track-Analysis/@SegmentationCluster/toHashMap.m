function hm = toHashMap(sc, varargin)
%converts segmentation cluster or group of segmentation clusters to 
%a hashmap or array of hashmaps that is friendly to yaml processing
%
%function hm = toHashMap(sc, varargin)

if (length(sc) > 1)
    hm = javaArray('java.util.HashMap', length(sc));
    for j = 1:length(sc)
        hm(j) = toHashMap(sc(j));
    end
else    
    no = convertSc(sc);
    hm = toHashMap(no);
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
if (isempty(kp) || isempty(kp(1).track) || isempty(kp(1).inds))
    np.filename = '';
    np.locInFile = [];
    newobj.knownPoints = np;
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

newobj.knownPoints = np;
