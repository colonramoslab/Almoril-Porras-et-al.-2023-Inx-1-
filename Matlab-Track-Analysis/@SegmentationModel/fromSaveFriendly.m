function sm2 = fromSaveFriendly(sm, eset)

fields = fieldnames(sm);
fields = setdiff(fields, {'eset', 'segmentationClusters'});
for j = 1:length(sm)
    sm2(j) = SegmentationModel();
    for k = 1:length(fields)
        sm2(j).(fields{k}) = sm(j).(fields{k});
    end
    sm2.segmentationClusters = sm.segmentationClusters.fromSaveFriendly(eset);
end

sm2.eset = eset;