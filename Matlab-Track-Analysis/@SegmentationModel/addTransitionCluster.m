function addTransitionCluster (sm, basecluster, newclusterlocation, newclustername, varargin) 
%function addTransitionCluster (sm, basecluster, newclusterlocation, newclustername,  varargin) 

 
pre = 0;
post = 1;

if (ischar(basecluster))
    bn = basecluster;
    basecluster = find(strcmpi(basecluster, {sm.segmentationClusters.name}));
    if (isempty(basecluster))
        disp(['can''t find cluster ' bn]);
        return;
    end
end

switch(lower(newclusterlocation))
    case {'before', 'pre'}
        dir = pre;
    case {'after', 'post'}
        dir = post;
    otherwise
        disp ('new location must be ''before'' or ''after'' current cluster');
        return;
end

sc = sm.segmentationClusters;
sc.condenseKnownPoints();

%create a new segmentation cluster
newsc = SegmentationCluster();
f = setdiff(fieldnames(newsc),{'name', 'knownPoints'});
for j = 1:length(f)
    newsc.(f{j}) = sc(basecluster).(f{j});
end
newsc.name = newclustername;

%find all the indices that precede or follow the the base cluster

kp = sc(basecluster).knownPoints;
valid = true(size(kp));
for j = 1:length(kp)
    if dir == pre
        start = kp(j).inds(find(diff(kp(j).inds) > 1) + 1);
        inds = start - 1;
    else
        stop = kp(j).inds(diff(kp(j).inds) > 1);
        inds = stop + 1;
    end
    newkp(j).track = kp(j).track;
    newkp(j).inds = inds;
    if (isempty(inds))
        valid(j) = false;
    end
end
newkp = newkp(valid);

newsc.knownPoints = newkp;

%now make sure these indices aren't marked in any other cluster
for j = 1:length(sc)
    for k = 1:length(newkp)
        ki = find([sc(j).knownPoints.track] == newkp(k).track);
        if (~isempty(ki))
            sc(j).knownPoints(ki).inds = setdiff(sc(j).knownPoints(ki).inds, newkp(k).inds);
        end
    end
end


%insert the new cluster into sm.segmentationClusters
if (dir == pre)
    scnew = [sc(1:(basecluster - 1)) newsc sc(basecluster:end)];
    
    nc = basecluster;
    bc = nc + 1;
else
    scnew = [sc(1:basecluster) newsc sc((basecluster + 1):end)];
    bc = basecluster;
    nc = bc + 1;
end

%adjust the allowed transitions
atold = sm.allowedTransitions;
atnew = zeros(size(atold) + [1 1]);
if (dir == pre)
    ncinds = [1:(nc -1) (bc + 1):length(atnew)];
    ocinds = [1:(basecluster-1) (basecluster + 1):length(atold)];
    atnew(ncinds, ncinds) = atold(ocinds, ocinds);
    atnew(ncinds,nc) = atold(ocinds,basecluster);
    atnew(bc,ncinds) = atold(basecluster,ocinds);
    atnew(bc,bc) = 1;
    atnew(nc,[nc bc]) = 1;
else
    ncinds = [1:(bc -1) (nc + 1):length(atnew)];
    ocinds = [1:(basecluster-1) (basecluster + 1):length(atold)];
    atnew(ncinds, ncinds) = atold(ocinds, ocinds);
    atnew(nc,ncinds) = atold(basecluster,ocinds);
    atnew(ncinds,bc) = atold(ocinds,basecluster);
    atnew([bc nc], nc) = 1;
    atnew(bc,bc) = 1;
end

sm.segmentationClusters = scnew;
sm.allowedTransitions = atnew;