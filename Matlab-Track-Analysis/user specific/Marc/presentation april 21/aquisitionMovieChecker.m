outres = [768 1024]; %output resolution of the movie; not tested except at [768 1024]
fn = '\\labnas1\Share\Ashley Extracted\Checkerboard\pretty movie\07_1715_tracks.bin';
fstub = '\\labnas1\Share\Phototaxis\Data\2010\Apr\07\1715\1715_';
chnm = '\\labnas1\Share\Phototaxis\Data\2010\Apr\07\1715\1715_background\1715_background0.jpg';
avifn = 'checkerMovie.avi'; %output movie name
movieCodec = 'Cinepak'; % choose from 'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE' or 'None'
%type 
%makeAvi = true;
%at command line before running script to make the movie


%load & autosegment tracks
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

%change the analysis rectange to have the same w:h ratio as outres
cx = 0.5 * (ar(1) + ar(3));
cy = 0.5 * (ar(2) + ar(4));
width = ar(3)-ar(1);
height = ar(4)-ar(2);
ratio = outres(1)/outres(2);
height = min(height, width*ratio);
width = min(width, height/ratio);
ar = round([cx-width/2 cy-height/2 cx+width/2 cy+height/2]);

%add a second coordinate system to give larva position in pixels of new
%lower resolution output image

if (~exist('gq', 'var'))
    [x,y] = meshgrid(1:2592, 1:1944);
    size(x)
    x = imresize(x(ar(2):ar(4), ar(1):ar(3)), outres);
    y = imresize(y(ar(2):ar(4), ar(1):ar(3)), outres);
    gq = GlobalQuantity();
    gq.xField = 'sloc';
    gq.xData(:,:,1) = uint16(x-1);
    gq.xData(:,:,2) = uint16(y-1);
    gq.derivationMethod = @GlobalQuantity.twoDinterpolation;
    [x,y] = meshgrid(1:outres(2), 1:outres(1));
    gq.yData(:,:,1) = uint16(x-1);
    gq.yData(:,:,2) = uint16(y-1);
    gq.fieldname = 'imloc';
    eset.executeExperimentFunction('addGlobalQuantity', gq);
end

%load in raw images & calculate a background
inds = 1:1:480;
if (~exist('imstack','var'))
    imstack = zeros(outres(1), outres(2), length(inds), 'uint8');
    
    for j = 1:length(inds)
        im = imread([fstub num2str(inds(j)) '.jpg']);
        im2 = imresize(im(ar(2):ar(4), ar(1):ar(3)), outres);
        imstack(:,:,j) = im2;
    end
    bak = min(imstack,[],3);
end

%load in checkerboard image
if (~exist('chim', 'var'))
    fullchim = imread(chnm);
    chim = imresize(fullchim(ar(2):ar(4), ar(1):ar(3)), outres);
end

close all;
figure(1); clf(1); set(gcf, 'Color', [0 0 0]);

%find only the tracks that start at the beginning
track = eset.expt.track;
tinds = find([track.startFrame] <= inds(1));
iml = zeros(length(tinds),2,length(inds));

%make a stack of point locations that matches up with the low res images
for j = 1:length(tinds)
    etinds = track(tinds(j)).getDerivedQuantity('mapPtsToInterped');
    etinds = interp1(1:length(etinds), etinds, inds, 'nearest', NaN);
    iml(j,:,isfinite(etinds)) = track(tinds(j)).getDerivedQuantity('imloc', false, 'inds', etinds);
end
end0 = 50; %how many frames to show before overlaying tracks
end1 = size(imstack,3) - 30; %when to stop overlaying tracks and start zooming
c = 'rgbcmyw'; c = [c c c c c c c];
msize = 20;
existsAndDefault('makeAvi', false);

if makeAvi
     aviobj = avifile(avifn, 'Quality', 90, 'Compression', 'Cinepak');
else
    aviobj = [];
end
%find out location of image to make movie
imshow(imstack(:,:,1));
u = get(gca, 'Units');
set(gca, 'Units', 'Pixels')
avirect = get(gca, 'Position');
set(gca, 'Units', u);


%show checkerboard & label
imshow(chim);
h = text(outres(2)/2, outres(1)/2, 'Projected Checkerboard Pattern', 'Color', [0 0 1], 'FontSize', 40, ...
    'HorizontalAlignment', 'center');
for j = 1:30
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
    pause(0.05);
end
delete(h);

%show first raw image & label
imshow(imstack(:,:,1));
h = text(outres(2)/2, outres(1)/2, 'Infrared Camera Image', 'Color', [1 0 0], 'FontSize', 40, ...
    'HorizontalAlignment', 'center');
for j = 1:30
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
    pause(0.05);
end

%show background subtracted raw image + checkerboard & label
imshow(imstack(:,:,1)-bak+chim/2);
h = text(outres(2)/2, outres(1)/2, 'Overlay', 'Color', [1 0 1], 'FontSize', 40, ...
    'HorizontalAlignment', 'center');
for j = 1:30
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
    pause(0.05);
end
delete(h)


%just larva on checkerboard
for j = 1:end0        
    imshow(imstack(:,:,j)-bak+chim/2);
    
    pause(0.05);
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
end

%larva on checkerboard with tracks overlayed
for j = end0:end1        
    imshow(imstack(:,:,j)-bak+chim/2); hold on
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
%plot all tracks to help me pick which ones to select; not needed now
%{
 imshow(imstack(:,:,j)-bak+chim/2); hold on
   
for k = 1:size(iml,1)
        plot (iml(k,1,j), iml(k,2,j), [c(k) 'o'], 'MarkerSize', msize)
        track(tinds(k)).plotPath('imloc', [c(k) '.'], 'inds', inds(end1):(inds(end1) + 240), 'MarkerSize', 2)
end
hold off
%}

%%
%select 4 tracks to zoom on, starting in upper left and going counter
%clockwise
si = [9 1 5 6];
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

%zoom in on specific animals
for j = ci:size(imstack,3)
    r = (j - ci) / (size(imstack,3) - ci);
    im4 = imstack(:,:,j) - bak + chim/2;
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
%%
%make 4 sub axes
inds2 = max(inds) + (1:120);

axloc = get(gca, 'position');

w = axloc(3)/2;
h = axloc(4)/2;
xl = axloc(1) + [0 0 w w];
yl = axloc(2) + [h 0 0 h];

clf(gcf);
for j = 1:4
    ax(j) = axes('Position', [xl(j) yl(j) w*.99 h*.99], 'Color', [0 0 0], 'XColor', tc(j), 'YColor', tc(j));
end
eset.expt.openDataFile;

for j = 1:4
    etinds = track2(j).getDerivedQuantity('mapPtsToInterped');
    interpinds{j} = interp1(etinds, 1:max(inds2), 'nearest', NaN);
end
imData.x = (1:size(fullchim,2)) - 1;
imData.y = (1:size(fullchim,1)) - 1;
imData.im = double(fullchim);


%scale = 0.99 * outres(2)/2 / 40;
scale = 1;
options = {'pretty', true, 'scale', scale, 'cropSize', [40 30],...
            'drawContour', true, 'contourColor', 'm-', 'LineWidth', 1, 'mhWidth', 3.5, 'underlayImData', imData};

%show extraction overlay on all 4 simultaneously        
figure(1);
for j = inds2    
    for k = 1:4        
        track2(k).pt(j).drawTrackImage([], 'Axes', ax(k), 'fid', eset.expt.fid , options{:});
        axis(ax(k),'ij');
        hold (ax(k), 'on')
        set(ax(k), 'CLim', [0 255])
        track2(k).plotPath('sloc', [tc(k) '.-'], 'Axes', ax(k), 'inds', interpinds{k}((j-50):1:j));
        hold (ax(k), 'off')
        set(ax(k), 'XTick', [], 'YTick', [],'Color', [0 0 0], 'Layer', 'top', 'Box', 'on', ...
            'XColor', tc(k), 'YColor', tc(k),'LineWidth',3,'YDir', 'reverse');
    end
    pause(0.05);
    if (~isempty(aviobj))
        F = getframe(gcf, avirect);
        aviobj = addframe(aviobj, F);
    end
end

%%
%show a signle track with speed, position, and body angle 
inds3 = max(inds2) + (1:480);
for j = 1:4
    set(ax(j), 'OuterPosition', [xl(j) yl(j) w*.99 h*.99], 'Color', [0 0 0]);
end
titleoptions = {'Color', 'w', 'FontSize', 14};
dpo = {'Color', 'w', 'LineWidth', 3};
dao = {'Color', [0 0 0], 'Layer', 'bottom', 'XTick', [], 'YTick', [], 'XColor', 'w', 'YColor', 'w', 'Box','on', 'LineWidth', 3};
lpo = {'Color', 'b', 'LineWidth', 2, 'imData', imData, ...
    'AxesOptions', {'Color', [0 0 0], 'Layer', 'bottom', 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k', 'Box','off'}};

figure(1);
aviobj = presentationMovie(track2(1), 'inds', inds3, 'AxesList', ax, options{:}, 'TitleOptions', titleoptions ,...
    'DataPlotOptions', dpo, 'DataAxesOptions', dao, 'ptbuffer', 200, 'LocPlotOptions', lpo,...
    'aviobj', aviobj, 'avirect', avirect, 'underlayImData', imData);
%%
%show individual panes with runs, reorientations, accepted headswings, and
%rejected headswings
figure(1);
aviobj = parallelMovies(track2, options{:}, 'AxesList', ax, 'aviobj', aviobj, 'avirect', avirect, 'npts', 240);

%%
if (~isempty(aviobj))
    aviobj = close(aviobj);
end
makeAvi = false;
