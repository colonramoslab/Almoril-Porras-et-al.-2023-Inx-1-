classdef TemperatureMaggotTrack < MaggotTrack
    %Temperature Maggot Track adds support for a global temperature
    %that changes with time; changes segmentation slightly to use
    %an adjusted speed (which you must specify using the tempToSA field)
    %rather than an absolute speed, since speed is a function of
    %temperature
    
    properties
        %tempToSA is a LUT that gives speed adjustment in terms
        %of temperature
        %SA = interp1(tempToSA(1,:), tempToSA(2,:), temp)
        %adjusted speed = SA*speed
        tempToSA = [];
    end
    
    methods %constructor
        function mt = TemperatureMaggotTrack (varargin)
            mt = mt@MaggotTrack(varargin);
            if ((nargin >= 1) && isa(varargin{1}, 'MaggotTrack'))
                mt.clone(varargin{1});
            end
            mt.so.speed_field = 'adjusted_speed';
        end
    end
    
    methods
        calculateDerivedQuantity (track, quantityNames, recalculate);
        %segmentTrack (track, maggotSegmentOptions);
    end
    methods (Static)
         varargout = validDQName (varargin)
    end
    
end

