fstub = '\\labnas1\Share\Phototaxis\Data\2010\Apr\01\1730\foreground';
start = 0:1440:5760;
stop = start + 1440;

for j = 1:length(start)
    inds = start(j):4:stop(j);
    tic
    [f,b] = makeForeground(fstub, inds);
    toc
    %imagesc(fg); colormap gray; colorbar vert
    fg(:,:,j) = f;
    bak(:,:,j) = b;
end

save ('\\labnas1\Share\Phototaxis\Data\2010\Apr\01\1730\marcfgextracted.mat');
%%
fg2 = fg;
fg2b = fg;
for j = 1:length(start)
    fg2(:,:,j) = imdilate(fg2(:,:,j), ones(100));
    fg2b(:,:,j) = conv2(gaussKernel(60), gaussKernel(60), fg2(:,:,j), 'same');
end
%%
foreground = max(fg2,[],3);
foreblur = max(fg2b,[],3);

foreblur = conv2(gausskernel(60), gausskernel(60), foreblur, 'same');


save ('\\labnas1\Share\Phototaxis\Data\2010\Apr\01\1730\marcfgextracted.mat');
%%
yinds = 1:1944;
xinds = 280:2280;

foreclip = max(fg,[],3);
foreclip = imdilate(foreclip, ones(150));
foreclip = double(foreclip(yinds, xinds));
inds = find(foreclip > 100);
[x,y] = meshgrid(xinds, yinds);
ydata = foreclip(inds);
order = 2;
lp = (order+1)*(order+2)/2;
v = myvander(x(inds), y(inds), lp);
p = pinv(v)*(ydata-mean(ydata));
v = myvander(x(:), y(:), lp);

fval = reshape(v*p+mean(ydata),size(foreclip));
figure(2);
imagesc(fval); colorbar('vert');