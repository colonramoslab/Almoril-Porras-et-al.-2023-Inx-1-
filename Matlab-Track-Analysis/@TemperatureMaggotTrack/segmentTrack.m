function segmentTrack (track, maggotSegmentOptions)
%function segmentTrack (track, maggotSegmentOptions)
%
if (nargin > 1)
    segmentTrack@MaggotTrack(track, maggotSegmentOptions);
else
    segmentTrack@MaggotTrack(track);
end