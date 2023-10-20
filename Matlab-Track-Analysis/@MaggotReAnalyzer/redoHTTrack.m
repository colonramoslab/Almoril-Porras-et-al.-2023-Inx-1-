function redoHTTrack(mra, track, varargin)
%function redoHTTrack(mra, track, varargin)
%

if (length(track) > 1)
    ts1 = tic;
    for j = 1:length(track)
        redoHTTrack(mra, track(j), varargin{:});
        disp ([num2str(j) ' - ' num2str(toc(ts1)) ' s']);
    end
    return;
end
redospine = true;
varargin = assignApplicable(varargin);

pt = track.pt;
nummidpts = 11;
ts = tic;
 lasttime = toc(ts);
for j = 1:length(pt)
    if (length(pt(j).contour) < 10)
        continue;
    end
    cpts = interp1(pt(j).contour',linspace(1, length(pt(j).contour), length(pt(j).contour)*mra.contourScale))';
    inds = findPointyEnds(cpts);
    d1 = sum((pt(j).head-cpts(:,inds(1))).^2) + sum((pt(j).tail-cpts(:,inds(2))).^2);
    d2 = sum((pt(j).head-cpts(:,inds(2))).^2) + sum((pt(j).tail-cpts(:,inds(1))).^2);
    if (d1 < d2)
        inds = inds([2 1]);
    end
    pt(j).head = cpts(:,inds(2));
    pt(j).tail = cpts(:,inds(1));
   
    if (redospine)
        if (toc(ts) - lasttime > 60)
            disp ([num2str(j) ' - ' num2str(toc(ts))]);
            lasttime = toc(ts);
        end
        ml = splitOutline(cpts, inds(1), inds(2));
        %newmid = snakeToFindMidline(pt(j).contour, pt(j).tail, pt(j).head);
%        [ml, c1, c2] = splitOutline(cpts, inds(1), inds(2));
 %       ml = refineMidline(ml, c1, c2);
        oldml = ml;
        ml = lowpass1D(ml, 0.5*(length(ml)/nummidpts),'padType','linear');
        ml(:,1) = oldml(:,1);
        ml(:,end) = oldml(:,end);
        newmid = resampleContour(ml, 'closed', false);
        pt(j).spine = interp1(newmid', linspace(1,length(newmid),nummidpts))';
        pt(j).mid = interp1(newmid', (length(newmid) + 1)/2)';
        pt(j).head = newmid(:,end);
        pt(j).tail = newmid(:,1);
    end
end
track.pt = pt;
    

