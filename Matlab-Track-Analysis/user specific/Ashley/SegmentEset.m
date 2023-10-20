function [eset] = SegmentEset(eset)
%function [eset] = SegmentEset(eset)
%
%after loading the eset using Load Multiple Files you might want to change
%the segmentation of the tracks.

for n=1:length(eset)
    expt=eset.expt(n);
    disp(['Analyzing file #' num2str(n)]);
    expt.executeTrackFunction('setSegmentSpeeds');
    expt.executeTrackFunction('segmentTrack');
    expt.executeTrackFunction('calculateDerivedQuantity', {'smid'})
    expt.executeTrackFunction('calculateDerivedQuantity', {'vmid'})
    eset.expt(n)=expt;
end

end

