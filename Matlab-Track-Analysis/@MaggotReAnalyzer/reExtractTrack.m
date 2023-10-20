function  reExtractTrack(mra, track, varargin)
%function  reExtractTrack(mra, track, varargin)

existsAndDefault('track', mra.track);
resegment = true;
varargin = assignApplicable(varargin);
pt2 = mra.rethreshold (track.pt(1), varargin{:});
pt2 = repmat(pt2, size(track.pt));
s = warning('off', 'all');
for j = 1:length(track.pt)
    pt2(j) = mra.rethreshold(track.pt(j), varargin{:});
    pt2(j) = mra.findHT (pt2(j), varargin{:});
end
warning(s);
track.pt = pt2;
mra.alignHTTrack(track);
track.recalculateDerivedQuantities;
if (resegment)
    track.setSegmentSpeeds();
    track.segmentTrack();
end
