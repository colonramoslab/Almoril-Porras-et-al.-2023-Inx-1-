classdef MartaTrack < MaggotTrack
    %Extension of MaggotTrack for Marta Zlatic data
    %main difference is how to load a file
    
    properties
        fname;
    end
    methods
        fixHTOrientation(track, varargin);
        markHTInvalid(track, thresh, varargin);
        valid = removeCollisionPoints(track, maxArea, varargin);        
        calculateDerivedQuantity(track, quantityNames, recalculate);
    end
    methods %constructor
        function mt = MartaTrack (varargin)
            mt = mt@MaggotTrack(varargin);
            mt.dr.interpTime = 0.025; %video taken at ~40 fps
            mt.dr.smoothTime = 0.05; %try less smoothing & see
            mt.dr.derivTime = 0.025; %1 pt derivative
            mt.so = MaggotSegmentOptions();
            mt.so.speed_field = 'smoothSpeed';
            if ((nargin >= 1) && isa(varargin{1}, 'MaggotTrack'))
                mt.clone(varargin{1});
            end
        end
    end %constructor
    methods (Static)
        mt = fromFile(fname);
        varargout = validDQName (varargin);
    end
    
end

