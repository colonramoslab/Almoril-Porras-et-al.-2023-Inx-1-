clear all
close all
clear classes
d = dir ('ispc.mat');
if (isempty(d))
    marcMacSetup
    ispc = true;
else
    marcPCSetup
    ispc = false;
end
fn = [basedir 'n2_g15_tracks.bin'];
timfn = [basedir 'n2_g15_.tim'];
tmpfn = [basedir 'n2_g15_.tmp'];
fstub = [basedir 'n2_g15_'];
tic; af = Experiment.fromFile (fn); toc
tic; af.addtime(timfn); toc
minFrames = 10;
af.track = af.track([af.track.npts] > minFrames);
tic; af.stitchTracks(10, 10); toc

[blah, I] = max([af.track.npts]);

track2 = af.reloadTrack(af.track(I));
return

im = imread([fstub '1000.jpg']);
figure(1);clf(1);imagesc(im); hold on; colormap gray

for j = 1:length(af.track)
    if (af.track(j).npts > 100)
         af.track(j).plotPath(); hold on
    else
         af.track(j).plotPath('r-'); hold on
    end
end