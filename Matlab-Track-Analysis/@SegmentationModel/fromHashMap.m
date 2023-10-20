function sm = fromHashMap(hm, varargin)
%function sm = fromHashMap(hm, varargin)

if isa(hm, 'java.util.ArrayList')
    it = hm.iterator;
    for j = 1:hm.size()
        try
            sm(j) = fromHashMap(it.next, varargin{:}); %matlab autoconverts basic types
        catch me
            disp('error');
            me.getReport
            break;
        end
    end
    return;
end

if ~isa(hm, 'java.util.HashMap')
    warning ('expecting a hash map or java arraylist of hashmaps');
    sm = [];
    return;
end
if (isa(varargin{1}, 'ExperimentSet'))
    eset = varargin{1};
    varargin = varargin(2:end);
else
    eset = [];
end
varargin = assignApplicable(varargin);
if (isempty(eset))
    error ('need an experiment set, so I can assign tracks');
end
sm = SegmentationModel();
sm.eset = eset;
sm.segmentationClusters = SegmentationCluster.fromHashMap(hm.get('segmentationClusters'), 'eset', eset, varargin{:});
hm.remove('segmentationClusters');
no = objFromHashmap(hm);
flist = intersect(fieldnames(no), fieldnames(sm));
for j = 1:length(flist)
    sm.(flist{j}) = no.(flist{j});
end
sm.allowedTransitions = reshape(sm.allowedTransitions, length(sm.segmentationClusters), []);
