function drawTrackImage(tp, camcalinfo, varargin)
%function drawTrackImage(tp, camcalinfo, varargin)
%
%draws a marta track point, spine, head tail
Axes = [];
varargin = assignApplicable(varargin);
if (isempty(Axes))
    Axes = gca;
end

if (all(isfinite(tp.contour)))
    plot (Axes, tp.spine(1,:), tp.spine(2,:), 'k-', tp.spine(1,:), tp.spine(2,:), 'y.', ...
        tp.head(1), tp.head(2), 'go', tp.tail(1), tp.tail(2), 'rh', tp.loc(1), tp.loc(2), 'kx',...
        tp.contour(1,[1:end 1]), tp.contour(2,[1:end 1]),'r-',...
        'MarkerSize', 10, 'LineWidth', 2);
else
    plot (Axes, tp.spine(1,:), tp.spine(2,:), 'k-', tp.spine(1,:), tp.spine(2,:), 'y.', ...
        tp.head(1), tp.head(2), 'go', tp.tail(1), tp.tail(2), 'rh', tp.loc(1), tp.loc(2), 'kx',...
        'MarkerSize', 10, 'LineWidth', 2);
end
