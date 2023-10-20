function regularizeSpine (track)
%function regularizeSpine (track)

pt = track.pt;
for j = 1:length(pt)
    pt(j).spine = resampleContour(pt(j).spine, 'closed', false);
end
track.pt = pt;

