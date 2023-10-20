classdef MartaTrackPoint < MaggotTrackPoint
    %an extension of Maggot Track Point to deal with data imported
    %from Marta Zlatic.
    %extra field is "spine" which represents an 11 point line down the
    %midpoint of the animal -- spine moved to MaggotTrackPoint 9/15/10
    
    properties
        
    end
    
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
        drawTrackImage(tp, camcalinfo, varargin);
        drawContourAndHead(pt, varargin);
        mtp = calculateCovariance(mtp);
    end
    
    methods
        function tp = MartaTrackPoint(varargin)
            tp = tp@MaggotTrackPoint(varargin{:});
            if (nargin >= 1) && (isa(varargin{1}, 'MartaTrackPoint'))
                op = varargin{1};
                flist = fieldnames(tp);
                for j = 1:length(flist)
                    tp.(flist{j}) = op.(flist{j});
                end
            end
        end%MaggotTrackPoint
   
    end
    
end

