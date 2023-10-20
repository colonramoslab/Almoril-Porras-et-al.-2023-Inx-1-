function sc = fromJavaArr(obj, eset)
%converts a java object (e.g. from sc.toJavaArr) into an array of
%segmentation clusters
%function arr = fromJavaArr(obj, eset)

if (~exist('eset', 'var'))
    disp ('eset is required to match points to tracks');
end

newobj = objFromHashmap(obj);

for j = 1:length(newobj)
    sc(j) = convertNewObj(newobj(j), eset);
end

%reverses this process:
%converts segmentation cluster to text friendly version by
%indentifying known points by file name and location in file, rather
%than by handle and index, which can change
function newobj = convertNewObj(newobj, eset)

sc = SegmentationCluster();
f = intersect(fieldnames(o2), fieldnames(sc));
f = setdiff(f, {'knownPoints', 'clustCov'});

for j = 1:length(f)
    sc.(f{j}) = newobj.(f{j});
end
newobj.clustCov = sc.clustCov(:);

sc.knownPoints = lookUpPoints(newobj.knownPoints, eset);


function newkp = lookUpPoints(kp, eset)



function ind = matchFileNames(fn, fnlist)
if any(sctrcmp(fn, fnlist))
    ind = find(sctrcmp(fn, fnlist));
    return;
end
if any(sctrcmpi(fn, fnlist))
    ind = find(sctrcmpi(fn, fnlist));
    return;
end
for j = 1:length(fnlist)
    [~,f] = fileparts(fnlist{j});
    short{j} = f;
end
[~,fn] = fileparts(fn);
ind = find(strcmpi(fn, short));
