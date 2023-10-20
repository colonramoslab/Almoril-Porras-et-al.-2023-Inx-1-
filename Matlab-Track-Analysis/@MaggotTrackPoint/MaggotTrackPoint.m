classdef MaggotTrackPoint < ImTrackPoint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        targetArea = -1;
        threshold = -1;
        htValid = false;
        head = [NaN NaN];
        mid = [NaN NaN];
        tail = [NaN NaN];
        contour = [NaN NaN];
        spine = [NaN NaN];
    end
    
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
        drawTrackImage(tp, camcalinfo, varargin);
        drawContourAndHead(pt, varargin);
    end
    
    methods
        function tp = MaggotTrackPoint(varargin)
            tp = tp@ImTrackPoint(varargin{:});
            if (nargin >= 1) && (isa(varargin{1}, 'MaggotTrackPoint'))
                op = varargin{1};
                flist = fieldnames(tp);
                for j = 1:length(flist)
                    tp.(flist{j}) = op.(flist{j});
                end
            end
        end%MaggotTrackPoint
   
    end
    %}
end

