outres = [768 1024];

fstub = '\\labnas1\Share\Phototaxis\Data\2010\Apr\02\1305\1305_';
fn = '\\labnas1\Share\Ashley Extracted\Temporal\ON OFF\starvation\02_1305_tracks.bin';
avifn = 'temporalMovie.avi';
 
if (~exist('eset', 'var'))
    eset = ExperimentSet.fromFiles(fn);
end

existsAndDefault('stitch',true);
if (stitch)
    fd = 3;
    md = 10;
    eset.expt(1).stitchTracks(fd, md, 'interactive', false);
    stitch = false;
end

existsAndDefault('segment', true);
if (segment)
    eset.executeTrackFunction('fixHTOrientation');
    eset.executeTrackFunction('setSegmentSpeeds');
    eset.executeTrackFunction('segmentTrack');
    segment = false;
end


ar = [300 100 2250 1900]; %analysis rectangle for analysis

%resize to fit ratio on output
cx = 0.5 * (ar(1) + ar(3));
cy = 0.5 * (ar(2) + ar(4));
width = ar(3)-ar(1);
height = ar(4)-ar(2);
ratio = outres(1)/outres(2);
height = min(height, width*ratio);
width = min(width, height/ratio);
ar = round([cx-width/2 cy-height/2 cx+width/2 cy+height/2]);

if (~exist('gq', 'var'))
    [x,y] = meshgrid(1:2592, 1:1944);
    size(x)
    x = imresize(x(ar(2):ar(4), ar(1):ar(3)), outres);
    y = imresize(y(ar(2):ar(4), ar(1):ar(3)), outres);
    gq = GlobalQuantity();
    gq.xField = 'sloc';
    gq.xData(:,:,1) = uint16(x);
    gq.xData(:,:,2) = uint16(y);
    gq.derivationMethod = @GlobalQuantity.twoDinterpolation;
    [x,y] = meshgrid(1:outres(2), 1:outres(1));
    gq.yData(:,:,1) = uint16(x);
    gq.yData(:,:,2) = uint16(y);
    gq.fieldname = 'imloc';
    eset.executeExperimentFunction('addGlobalQuantity', gq);
end

inds = 1:1:240;
if (~exist('imstack','var'))
    imstack = zeros(outres(1), outres(2), length(inds), 'uint8');

    for j = 1:length(inds)
        im = imread([fstub num2str(inds(j)) '.jpg']);
        im2 = imresize(im(ar(2):ar(4), ar(1):ar(3)), outres);
        imstack(:,:,j) = im2;
    end
end


close all;
figure(1); clf(1); set(gcf, 'Color', [0 0 0]);

track = eset.expt.track;
tinds = find([track.startFrame] <= inds(1));
iml = zeros(length(tinds),2,length(inds));
for j = 1:length(tinds)
    etinds = track(tinds(j)).getDerivedQuantity('mapPtsToInterped');
    etinds = interp1(1:length(etinds), etinds, inds, 'nearest', NaN);
%    etinds = etinds(inds);
%    size(track(tinds(j)).getDerivedQuantity('imloc', false, 'inds', etinds))
    iml(j,:,isfinite(etinds)) = track(tinds(j)).getDerivedQuantity('imloc', false, 'inds', etinds);
end
end0 = 50;
end1 = size(imstack,3) - 30;
c = 'rgbcmyw'; c = [c c c c c c c];
msize = 20;
existsAndDefault('makeAvi', false);
if makeAvi
     aviobj = avifile(avifn, 'Quality', 90, 'Compression', 'Cinepak');
     
else
    aviobj = [];
end
imshow(imstack(:,:,1));
u = get(gca, 'Units');
set(gca, 'Units', 'Pixels')
avirect = get(gca, 'Position');
set(gca, 'Units', u);
for j = 1:end0        
    imshow(imstack(:,:,j)); 
    
    pause(0.05);
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
end

for j = end0:end1        
    imshow(imstack(:,:,j)); hold on
    for k = 1:size(iml,1)
        plot (iml(k,1,j), iml(k,2,j), [c(k) 'o'], 'MarkerSize', msize)
        track(tinds(k)).plotPath('imloc', [c(k) '.'], 'inds', inds(1:j), 'MarkerSize', 2)
    end
    hold off
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
    pause(0.05);
end

%%
%select 4 tracks to zoom on, starting in upper left and going counter
%clockwise
si = [6 4 1 9];
track2 = track(tinds(si));
tc = c(si);
start = iml(si,:,end1);
xx = outres(2)/4;
yy = outres(1)/4;
finish = [xx yy; xx 3*yy; 3*xx 3*yy; 3*xx yy];
startsize = [40 30];
endsize = [outres(2) outres(1)]/2 * .99;

figure(1); clf(1);
slide = @(x1,x2,r) x1 + (x2 - x1)*r;
ci = end1;
numleft = size(imstack,3) - ci;
for j = ci:size(imstack,3)
    r = (j - ci) / (size(imstack,3) - ci);
    im4 = imstack(:,:,j);
    for k = 1:length(si)
        [im4, z] = drawZoomBox(imstack(:,:,j), im4, iml(si(k),:,j), startsize, ...
            slide(iml(si(k),:,j), finish(k,:), r), slide(startsize, endsize, r));
        zb{k} = z;
    end
    imshow (im4);
    hold on
    for k = 1:length(zb)
        plot (zb{k}{:}, [tc(k) '-'], 'LineWidth', 3);
    end
    hold off
    pause(0.05);
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
end

inds2 = max(inds) + (1:120);

axloc = get(gca, 'position');

w = axloc(3)/2;
h = axloc(4)/2;
xl = axloc(1) + [0 0 w w];
yl = axloc(2) + [h 0 0 h];
%%
clf(gcf);
for j = 1:4
    ax(j) = axes('Position', [xl(j) yl(j) w*.99 h*.99], 'Color', [0 0 0], 'XColor', tc(j), 'YColor', tc(j));
end
eset.expt.openDataFile;

for j = 1:4
    etinds = track2(j).getDerivedQuantity('mapPtsToInterped');
    interpinds{j} = interp1(etinds, 1:max(inds2), 'nearest', NaN);
end
%%

%scale = 0.99 * outres(2)/2 / 40;
scale = 1;
options = {'pretty', true, 'scale', scale, 'cropSize', [40 30],...
            'drawContour', true, 'contourColor', 'm-', 'LineWidth', 1, 'mhWidth', 3.5};
%%
for j = inds2    
    for k = 1:4        
        track2(k).pt(j).drawTrackImage([], 'Axes', ax(k), 'fid', eset.expt.fid , options{:});
        axis(ax(k),'ij');
        hold (ax(k), 'on')
        set(ax(k), 'CLim', [0 255])
        track2(k).plotPath('sloc', [tc(k) '.-'], 'Axes', ax(k), 'inds', interpinds{k}((j-50):1:j));
        hold (ax(k), 'off')
        set(ax(k), 'XTick', [], 'YTick', [],'Color', [0 0 0], 'Layer', 'top', 'Box', 'on', 'XColor', tc(k), 'YColor', tc(k),'LineWidth',3);
    end
    pause(0.05);
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
end

%%
inds3 = max(inds2) + (1:480);
for j = 1:4
    set(ax(j), 'OuterPosition', [xl(j) yl(j) w*.99 h*.99], 'Color', [0 0 0]);
end
titleoptions = {'Color', 'w', 'FontSize', 14};
dpo = {'Color', 'w', 'LineWidth', 3};
dao = {'Color', [0 0 0], 'Layer', 'bottom', 'XTick', [], 'YTick', [], 'XColor', 'w', 'YColor', 'w', 'Box','on', 'LineWidth', 3};
lpo = {'Color', 'w', 'LineWidth', 2, 'AxesOptions', {'Color', [0 0 0], 'Layer', 'bottom', 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k', 'Box','off'}};
aviobj = presentationMovie(track2(1), 'inds', inds3, 'AxesList', ax, options{:}, 'TitleOptions', titleoptions ,...
    'DataPlotOptions', dpo, 'DataAxesOptions', dao, 'ptbuffer', 200, 'LocPlotOptions', lpo,...
    'aviobj', aviobj, 'avirect', avirect);
%%
figure(1);
aviobj = parallelMovies(track2, options{:}, 'AxesList', ax, 'aviobj', aviobj, 'avirect', avirect, 'npts', 240);

%%
if (~isempty(aviobj))
    aviobj = close(aviobj);
end
makeavi = false;
