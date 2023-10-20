function drawContourAndHead(pt, varargin)
%function drawTrackImage (pt)
%@MaggotTrackPoint

Axes = [];
varargin = assignApplicable(varargin);
if (isempty(Axes))
    Axes = gca;
end

if (length(pt) > 1)
    ih = ishold(Axes);
    for j = 1:length(pt)
        pt(j).drawContourAndHead('Axes', Axes, varargin{:});
        hold(Axes, 'on');
    end
    if (~ih)
        hold(Axes, 'off');
    end
    return;
end
    
    

offset = [0;0];
contourColor = 'k-';

varargin = assignApplicable(varargin);

h = pt.head + offset;
m = pt.mid + offset;
t = pt.tail + offset;
c = pt.contour + repmat(offset, 1, length(pt.contour));
c(:,end+1) = c(:,1); 


plot (Axes, h(1),h(2),'g*', t(1),t(2),'rh', [t(1),m(1),h(1)], [t(2),m(2),h(2)], 'y-', c(1,:), c(2,:), contourColor, varargin{:});
