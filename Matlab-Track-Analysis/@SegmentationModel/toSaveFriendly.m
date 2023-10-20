function sm2 = toSaveFriendly(sm, varargin) 
%function sm2 = toSaveFriendly(sm, varargin) 

fields = fieldnames(sm);
fields = setdiff(fields, {'eset', 'segmentationClusters'});
for j = 1:length(sm)
    sm2(j) = SegmentationModel();
    for k = 1:length(fields)
        sm2(j).(fields{k}) = sm(j).(fields{k});
    end
    sm2.segmentationClusters = sm.segmentationClusters.toSaveFriendly();
end
