classdef ImTrackPoint < TrackPoint
    %Augments TrackPoint by including an excerpted image
    %   Detailed explanation goes here
    
    properties
        imOffset = int16([0 0]);
        imData = uint8([]);
    end
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
        drawTrackImage(tp, camcalinfo, varargin)
    end
    
    methods %constructor
        function tp = ImTrackPoint(varargin)
            tp = tp@TrackPoint(varargin{:});            
            if (nargin >= 1) && (isa(varargin{1}, 'ImTrackPoint'))
                op = varargin{1};
                flist = fieldnames(tp);
                for j = 1:length(flist)
                     tp.(flist{j}) = op.(flist{j});
                end
            end
        end%trackpoint
    end
    
end

