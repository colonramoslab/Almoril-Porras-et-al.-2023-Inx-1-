classdef MaggotTrack < Track
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Transient = true, AbortSet = true)
        headSwing;
    end
    
    methods %constructor
        function mt = MaggotTrack (varargin)
            mt = mt@Track(varargin);
            mt.dr.smoothTime = 0.5; %maggots move faster & need less smoothing than worms
            mt.so = MaggotSegmentOptions();
            if ((nargin >= 1) && isa(varargin{1}, 'MaggotTrack'))
                mt.clone(varargin{1});
            end
        end
    end %constructor
    
    methods %access
        function set.headSwing(obj, value)
            if (~isempty(value) && isa(value, 'TrackPart'))
                [value.track] = deal(obj);
            end
            obj.headSwing = value;
        end
    end
    
    methods(Static)
        varargout = validDQName (varargin)
    end
    methods 
        recalculateDerivedQuantities(track, varargin); %recalculate all already derived quantities
        calculateDerivedQuantity(track, quantityName, recalculate);
        setSegmentSpeedsByPercentile (track, maggotSegmentOptions, stopPctl, startPctl);
        segmentTrack (track, maggotSegmentOptions);
        setSegmentSpeeds (track, mso);
        plotSegmentation (track, varargin);
        fixHTOrientation(track, varargin);
        prettyMovie(track, varargin);
    end
    
        
    
    
end

