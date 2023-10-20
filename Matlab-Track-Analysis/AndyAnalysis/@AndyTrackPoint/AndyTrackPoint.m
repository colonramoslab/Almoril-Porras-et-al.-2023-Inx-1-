classdef AndyTrackPoint < MartaTrackPoint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        mcdf = Mcd_Frame;
    end
    
    methods %constructor
        function atp = AndyTrackPoint(varargin)
            ismcdffun = @(m) isa(m, 'Mcd_Frame');
            ismcdf = cellfun(ismcdffun, varargin);
            arglist = varargin(~ismcdf);
            mcd = varargin(ismcdf);
            atp = atp@MartaTrackPoint(arglist{:});
            if (~isempty(mcd)) 
                atp.mcdf = mcd{1};
                atp = atp.populateMaggotFields;
            end
        end
    end
    methods
       atp = populateMaggotFields(atp);
    end
    methods (Static)
        function mmpp = mmPerPixel() 
            mmpp = 0.005;
        end
        function mmpsu = mmPerStageUnit()
            mmpsu = 0.0001;
        end
    end
    
end

