classdef CameraCalibration < handle
    %maps points from camera to real points
    
    properties
        realx;
        realy;
        camx;
        camy;
    end
    methods %constructor
        function cc = CameraCalibration(varargin)
            %cc = CameraCalibration(checkerim)
            %cc = CameraCalibration (realx, realy, camx, camy);
            switch(length(varargin))
                case 0
                    return;
                case 1
                    %assume it's an image of a checkerboard
                    [realx, realy, camx, camy] = calibrateCheckerboard(varargin{1});
                    cc = CameraCalibration(realx, realy, camx, camy);
                case 4
                    cc.realx = varargin{1};
                    cc.realy = varargin{2};
                    cc.camx = varargin{3};
                    cc.camy = varargin{4};
                    cc.setTSI();
                otherwise
                    try
                        [realx, realy, camx, camy] = calibrateCheckerboard(varargin{:});
                        cc = CameraCalibration(realx, realy, camx, camy);
                    catch me
                        disp(me.getReport);
                        disp('CameraCalibration(im) or CameraCalibration(realx, realy, camx, camy)');
                    end
                    return;
            end
        end    
    end
    
    methods
        [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin);
        camim = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin);
        campts = camPtsFromRealPts(cc, realpts);
        realpts = realPtsFromCamPts(cc, campts);
        magfactor = pixelsPerRealUnit (cc);
        magfactor = realUnitsPerPixel (cc);
        setTSI(cc); %set tri scattered interpolation
    end
    
    properties (SetAccess = protected)
        %tri scattered interpolants
        c2rX; 
        c2rY;
        r2cX;
        r2cY;
    end
end

