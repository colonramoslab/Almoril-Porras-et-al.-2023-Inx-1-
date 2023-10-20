classdef OldMaggotTrackPoint < MaggotTrackPoint
    % This just exists so we can load old maggot track points successfuly
    % (do not have midline saved in file)
    % do not construct normally -- fromFile returns a MaggotTrackPoint
    
    properties
       
    end
    
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
    end
    
    methods
        function tp = OldMaggotTrackPoint(varargin)
            tp = tp@MaggotTrackPoint(varargin{:});
        end%MaggotTrackPoint
   
    end
    %}
end

