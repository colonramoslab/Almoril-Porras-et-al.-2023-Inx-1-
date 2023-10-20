function [im, xaxis, yaxis] = trackTimeLapse (track, varargin)
%function [im, xaxis, yaxis] = trackTimeLapse (track, varargin)
%
%varargin: indsRange = [];
%underlayPath = false;
%indices = [];  if specified, use these points exactly to draw track
frameInterval = 100;
indsRange = [];
indices = [];
underlayPath = false;
varargin = assignApplicable(varargin);

sl = max(track.getDerivedQuantity('spineLength', false));
if (isempty(sl) || sl <= 0 || ~isfinite(sl))
    inds = 1:frameInterval:length(track.getDerivedQuantity('eti'));
end
pl = track.getDerivedQuantity('pathLength');

if (isempty(indsRange))
    if (~isempty(indices))
        indsRange = [min(indices) max(indices)];
    else
        indsRange = [1 length(pl)];
    end
end
start = pl(min(indsRange));
stop = pl(max(indsRange));
    
[pl, I] = unique(pl);
ds = sl*1.01;
inds = interp1(pl, I, start:ds:stop, 'nearest');

hs = track.headSwing;
if (~isempty(hs))
    mI = [hs.maxInd];
    mI = sort(mI(mI > indsRange(1) & mI < indsRange(end)));
else 
    mI = [];
end
if (length(mI) > 1)
    %space the headswings out by a distance greater than ds
    while (any(diff(pl(mI)) < ds))
        d = diff(pl(mI));
        [~,I] = min(d);
        if (I == 1)
            tokill = I + 1;
        else
            if (I == length(d))
                tokill = I;
            else
                if (d(I - 1) > d(I + 1))
                    tokill = I + 1;
                else
                    tokill = I;
                end
            end
        end
        mI = setdiff(mI, mI(I));
    end
    fixedpts = mI;
    if (pl(mI(1)) - pl(indsRange(1)) > ds)
        fixedpts = [indsRange(1) fixedpts];
    end
    if (pl(indsRange(end)) - pl(mI(end)) > ds)
        fixedpts = [fixedpts indsRange(end)];
    end
    inds = fixedpts(1);
    for j = 2:length(fixedpts)
        temppl = pl(fixedpts(j-1):fixedpts(j)) - pl(fixedpts(j-1));
        npts = floor(temppl(end) / ds);
        plinds = linspace(temppl(1), temppl(end), npts + 1);
        inds = [inds, interp1(temppl, fixedpts(j-1):fixedpts(j), plinds(2:end), 'nearest')];
    end
%{        
        if (pl(fixedpts(j)) - pl(fixedpts(j-1)) < 2*ds)
            inds = 
        
    dI = median(diff(inds)) * 0.2;
    for j = 1:length(inds)
        ii = find(mI + dI >= inds(j) & mI - dI < inds(j), 1, 'first');
        if (~isempty(ii))
            %inds(j)
            %mI(ii)
            inds(j) = mI(ii);
        end
    end
    %}
end

if (~isempty(indices))
    inds = indices;
end

inds = track.getDerivedQuantity('mapinterpedtopts',false,inds);



track.expt.openDataFile;
fid = track.expt.fid;
il = track.getDerivedQuantity('iloc');
ll = min(il(:,indsRange(1):indsRange(end)), [], 2); 
ur = max(il(:,indsRange(1):indsRange(end)), [], 2); 
if (~isempty(track.expt.camcalinfo))
    il = round(track.expt.camcalinfo.camPtsFromRealPts(il));
    ll = floor(track.expt.camcalinfo.camPtsFromRealPts(ll));
    ur = ceil(track.expt.camcalinfo.camPtsFromRealPts(ur));
end
pt = (track.pt(1));
fseek(fid, pt.locInFile, -1);
pt = pt.fromFile(fid, true, true, []);
ll = ll - max(size(pt.imData));
ur = ur + max(size(pt.imData));

xaxis = ll(1):ur(1);
yaxis = ll(2):ur(2);
im = double(zeros(ur(2)-ll(2) + 1, ur(1)-ll(1) + 1));


for j = inds
    pt = (track.pt(j));
    fseek(fid, pt.locInFile, -1);
    pt = pt.fromFile(fid, true, true, []);
    xinds = floor(double(pt.imOffset(1))-ll(1) + (1:size(pt.imData,2)));
    yinds = floor(double(pt.imOffset(2))-ll(2) + (1:size(pt.imData,1)));
    im(yinds, xinds) = max(im(yinds,xinds), double(pt.imData));
end    
size(im)
if (underlayPath)
    im = ((im - min(im(:)))/(max(im(:)) - min(im(:))));
    im = repmat(im, [1 1 3]);
    ind = sub2ind(size(im), il(2,indsRange(1):indsRange(end))-ll(2) + 1, il(1,indsRange(1):indsRange(end)) - ll(1) + 1, ones(size(il(1,indsRange(1):indsRange(end)))));
    im(ind) = 1;
end
size(im)