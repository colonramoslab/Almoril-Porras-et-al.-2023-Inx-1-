function sc2 = fromSaveFriendly(sc, eset)
%function sc2 = fromSaveFriendly(sc, eset)

for j = 1:length(sc)
    sc2(j) = convertNewObj(sc(j), eset);
end

%reverses this process:
%converts segmentation cluster to text friendly version by
%indentifying known points by file name and location in file, rather
%than by handle and index, which can change
function sc = convertNewObj(newobj, eset)

sc = SegmentationCluster();
f = intersect(fieldnames(newobj), fieldnames(sc));
f = setdiff(f, {'knownPoints'});

for j = 1:length(f)
    sc.(f{j}) = newobj.(f{j});
end

sc.knownPoints = lookUpPoints(newobj.knownPoints, eset);

function newkp = lookUpPoints(kp, eset)

if (isempty(eset) || isempty(kp) || isempty(kp(1).locInFile) || isempty(kp(1).filename))
    newkp = [];
    return;
end

for j = 1:length(kp)
    ind = matchFileNames(kp(j).filename, {eset.expt.fname});
    if isempty(ind)
        warning('sc:fhm', ['could not match file name: ' kp(j).filename]);
        kp(j)
        continue;
    end
    [t,inds] = matchLocToTrack(eset.expt(ind), kp(j).locInFile);
    for k = 1:length(t)
        newkp(k).track = t(k);
        newkp(k).inds = inds{k};
    end
end
if (~exist('newkp', 'var'))
    newkp = [];
end





function ind = matchFileNames(fn, fnlist)

if any(strcmp(fn, fnlist))
    ind = find(strcmp(fn, fnlist));
    return;
end
if any(strcmpi(fn, fnlist))
    ind = find(strcmpi(fn, fnlist));
    return;
end
for j = 1:length(fnlist)
    ftemp = fnlist{j};
    if (ispc)
        ftemp(ftemp == '/') = '\';
    else
        ftemp(ftemp == '\') = '/';
    end
    [~,f] = fileparts(ftemp);
    short{j} = f;
end
if (ispc)
    fn(fn == '/') = '\';
else
    fn(fn == '\') = '/';
end
    
[~,fn] = fileparts(fn);
ind = find(strcmpi(fn, short));

function [t,inds] = matchLocToTrack(expt, lif) 
    track = expt.track;
    pt = [track.pt];
    locInFile = [pt.locInFile];
    tnum = zeros(size(locInFile));
    ind0 = 0;
    for j = 1:length(track)
        inds = ind0 + (1:length(track(j).pt));
        tnum(inds) = j;
        ind0 = inds(end);
    end
    %[C,IA,IB] = INTERSECT(A,B) also returns index vectors IA and IB
    %such that C = A(IA) and C = B(IB).

    [lif, ~, I] = intersect(lif, locInFile);
    tnum = tnum(I);
    tlist = unique(tnum);
    clear inds;
    for j = 1:length(tlist)
        t(j) = track(tlist(j));
        pt = t(j).pt;
        [~,~,ptInds] = intersect(lif(tnum == tlist(j)), [pt.locInFile]);
        ii = ismember(t(j).getDerivedQuantity('mapinterpedtopts'), ptInds);
        inds{j} = find(ii);
    end
