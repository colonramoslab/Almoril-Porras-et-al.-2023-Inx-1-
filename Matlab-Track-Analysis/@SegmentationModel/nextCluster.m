function cnum = nextCluster(sm, currentcluster)
%
%most likely next cluster according to allowed transitions & transition
%matrix


if (ischar(currentcluster))
    bn = currentcluster;
    currentcluster = find(strcmpi(currentcluster, {sm.segmentationClusters.name}));
    if (isempty(currentcluster))
        disp(['can''t find cluster ' bn]);
        return;
    end
end

next = sm.allowedTransitions(currentcluster,:);
next(currentcluster) = 0;
next = find(next);
if (length(next) == 1)
    cnum = next;
    return;
end
if (any(size(sm.hmm_transmat) ~= length(sm.segmentationClusters)))
    [~,I] = max([sm.segmentationClusters(next).priorProbability]);
    cnum = next(I);
else
    [~,I] = max(sm.hmm_transmat(currentcluster, next));
    cnum = next(I);
end
