function tracks = getTracks(s)
%function tracks = getTracks(s)
%converts locations in simulation to tracks
pt = TrackPoint();
pt = repmat (pt, [1 s.npts]);
temp = num2cell(0:(s.npts - 1));
[pt.ind] = temp{:};
temp = num2cell(s.timestep * (0:(s.npts - 1)));
[pt.et] = temp{:};
for j = 1:s.ntracks
    temp = num2cell(squeeze(s.loc(j,:,:)), 1);
    [pt.loc] = temp{:};
    
    tracks(j) = Track();
    tracks(j).pt = pt;
    tracks(j).npts = s.npts;
    tracks(j).startFrame = 0;
    tracks(j).endFrame = s.npts - 1;
end