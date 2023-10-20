classdef MartaSegmentOptions < MaggotSegmentOptions
    %UNTITLED21 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hunch.formula = @(spineLength, spineWidth) (spineWidth - spineLength) > 4;
    end
    
    methods
        function mso = MartaSegmentOptions
               mso.speed_field = 'smoothSpeed'; %name of field that contains speed to use (normally 'speed')
               mso.stop_speed_cut = 0.7; %if speed < stop_speed_cut, end a run
               mso.start_speed_cut = 1.1; %if speed > start_speed_cut && vel_dp > aligned_dp, run can start
               mso.minRunTime = 1; %if run is less than this many seconds, discard 
        end
    end
    
end

