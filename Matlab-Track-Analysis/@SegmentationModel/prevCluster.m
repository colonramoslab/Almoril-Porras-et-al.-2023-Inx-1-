function cnum = prevCluster(sm, currentcluster)
%
%most likely prev cluster according to allowed transitions & transition
%matrix


if (ischar(currentcluster))
    bn = currentcluster;
    currentcluster = find(strcmpi(currentcluster, {sm.segmentationClusters.name}));
    if (isempty(currentcluster))
        disp(['can''t find cluster ' bn]);
        return;
    end
end

prev = sm.allowedTransitions(:,currentcluster);
prev(currentcluster) = 0;
prev = find(prev);
if (length(prev) == 1)
    cnum = prev;
    return;
end
if (any(size(sm.hmm_transmat) ~= length(sm.segmentationClusters)))
    [~,I] = max([sm.segmentationClusters(prev).priorProbability]);
    cnum = prev(I);
else
    [~,I] = max(sm.hmm_transmat(prev, currentcluster));
    cnum = prev(I);
end
