function playMovie(track, switchtime, varargin)
%@MaggotTrack
%playMovie(track, varargin)
%enter options as pairs, caps matter
%options, with defaults
%
%ptbuffer = 1000;
%delayTime = 0.05;
%axisSize (size of image or 120 * mean speed)
%inds = 1:length(track.pt);
%startLoc = []; > if startLoc & stopLoc are not empty, we run the movie
%between these two points
%stopLoc = []; >
%

darkmap = gray(256);
lightmap = darkmap(end:-1:1,:);

ptbuffer = 400;
delayTime = 0.05;
axisSize = max(size(track.pt(1).imData));
if (axisSize <= 0)
    axisSize = 100;
end
inds = 1:length(track.pt);
startLoc = [];
stopLoc = [];
pvpmod(varargin);
if (~isempty(startLoc) && ~isempty(stopLoc))
    [blah,s] = track.nearestPoint (startLoc);
    [blah,e] = track.nearestPoint (stopLoc);
    if (s > e)
        inds = e:s;
    else
        inds = s:e;
    end
end
loc = [track.pt.loc];
sloc = track.getDerivedQuantity('sloc');
sind = track.getDerivedQuantity('mapptstointerped');
track.calculateDerivedQuantity({'sbodytheta', 'speed', 'vel_dp'});
stheta = track.getDerivedQuantity('sbodytheta');
eti = track.getDerivedQuantity('eti');
head = track.getDerivedQuantity('shead');
so = track.so;

pt = [track.pt];
pauseCount = 1;
for j = inds
    ts = tic;
    hold off; cla()
    pt(j).drawTrackImage; hold on
    if (mod(floor(eti(sind(j))/switchtime), 2) == 0)
        colormap(lightmap);
        bc = [1 1 1];
    else
        colormap(darkmap);
        bc = [0 0 0];
    end
    set(gca, 'XTick', [], 'YTick', [], 'Color', bc);
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    plot (sloc(1,sstart:send), sloc(2,sstart:send), 'b.-');
    plot (sloc(1,sind(j)), sloc(2,sind(j)), 'b.', 'MarkerSize', 20);
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    axis equal; 
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
     if (~isempty(track.run))
        t = [];
        contourColor = 'm-';
        if (track.isrun(sind(j)))
            t = [t 'run '];
            contourColor = 'c-';
        end
        if (any([track.reorientation.inds] == sind(j)))
            t = [t 'reorientation '];
        end
        if (any([track.headSwing.inds] == sind(j)))
            t = [t 'headsweep '];
            I = find([track.headSwing.startInd] <= sind(j) & [track.headSwing.endInd] >= sind(j));
            if (~isempty(I))
                if (track.headSwing(I).accepted)
                    t = [t 'accepted '];
                    contourColor = 'g-';
                else
                    t = [t 'rejected '];
                    contourColor = 'r-';
                end
                ei = track.headSwing(I).endInd;
                si = track.headSwing(I).startInd;
                si = find(sign(1:si) ~= sign(stheta(si)) | abs(stheta(1:si)) < so.headswing_stop, 1, 'last');
                if (isempty(si))
                    si = track.headSwing(I).startInd;
                end
                plot (head(1,si:ei), head(2,si:ei), contourColor, 'LineWidth', 2);
            end
        end
        c = pt(j).contour;
        c(:,end+1) = c(:,1);
        plot (c(1,:), c(2,:), contourColor, 'LineWidth', 2);
       % t = [t num2str(eti(sind(j))) ' s elapsed'];
        title (t); embiggen()
     end
    hold off
    set(gca, 'XTick', [], 'YTick', []);
    timeleft = delayTime - toc(ts);
    if (timeleft > 0)
        pause(timeleft);
    else
        refresh(gcf);
    end
        
end



