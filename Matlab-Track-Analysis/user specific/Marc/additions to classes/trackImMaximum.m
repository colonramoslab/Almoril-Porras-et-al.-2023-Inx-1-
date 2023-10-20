function im = trackImMaximum (track, im, lowerLeft, varargin)
%function im = trackImMaximum (track, im, lowerLeft, varargin)

existsAndDefault('im', zeros([1944 2592]));
existsAndDefault('lowerLeft', [0 0]);
if (length(track) > 1)
    updateEvery = 1;
    varargin = assignApplicable(varargin);
    
    ts1 = tic;
    for j = 1:length(track)
        im = trackImMaximum(track(j), im, lowerLeft, varargin{:});
        if (updateEvery > 0 && mod(j, updateEvery) == 0)
            imagesc(im);
            title ([num2str(j) '/' num2str(length(track)) ' -- ' num2str(toc(ts1)) ' s']);
            
            pause(0.01);
        end
    end
    return
end

frameInterval = 10;
displacement = false;
varargin = assignApplicable(varargin);

inds = 1:frameInterval:length(track.getDerivedQuantity('eti'));

inds = track.getDerivedQuantity('mapinterpedtopts',false,inds);



track.expt.openDataFile;
fid = track.expt.fid;

ll = lowerLeft;
ptarr = track.pt(inds);
if (displacement)
    pt = ptarr(1);
    fseek(fid, pt.locInFile, -1);
    pt = pt.fromFile(fid, true, true, []);
    offsetX = round(size(im, 2)/2 - pt.loc(1));
    offsetY = round(size(im, 1)/2 - pt.loc(2));
else
    offsetX = 0;
    offsetY = 0;
end
for j = 1:length(ptarr)
    pt = ptarr(j);
    fseek(fid, pt.locInFile, -1);
    pt = pt.fromFile(fid, true, true, []);
    xinds = floor(double(pt.imOffset(1))-ll(1) + (1:size(pt.imData,2))) + offsetX;
    yinds = floor(double(pt.imOffset(2))-ll(2) + (1:size(pt.imData,1))) + offsetY;
    validindsX = xinds >= 1 & xinds <= size(im,2);
    validindsY =  yinds > 1 & yinds <= size(im,1);
    xinds = xinds(validindsX);
    yinds = yinds(validindsY);
    im(yinds, xinds) = max(im(yinds,xinds), double(pt.imData(validindsY, validindsX)));
end    
