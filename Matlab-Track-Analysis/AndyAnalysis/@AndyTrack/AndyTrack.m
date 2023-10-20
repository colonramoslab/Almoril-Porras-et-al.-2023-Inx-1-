classdef AndyTrack < MartaTrack
    %Track to analyze andy's dlp data sets
    %I extended marta track instead of maggot track, because I think the
    %faster temporal frequency in those experiments is likely reflected
    %here
    
    properties
    end
    
    methods (Static)
        at = fromFile(fname);
    end
    
end

