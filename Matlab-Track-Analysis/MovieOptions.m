classdef MovieOptions 
    
    properties
     makeAvi = false; %whether to make an avi from this movie
     aviName = 'none'; %name of avi file
     avioptions = {}; %options to pass to avifile when creating avi
                        %see help avifile for list
                        %example  avioptions = {'FPS', 10} to set frame
                        %rate to 10
     closeAviWhenDone = true; %if true, close the movie if you made it
     addToAvi = false; %if true, add to an existing avi file, aviName, avioptions not used 
     avi = [];        
    
     aviResolution = [1024 768];
     jpegResolution = [1024 768];
     jpegstub = '';
     makeImStack = false;
     imstackIndex = 0;
     
     
     
     backgroundColor = [0 0 0];
   
     figHandle = [];
     Axes = []; %axes into which to play movie
     datafields = {'speed', 'deltatheta', 'covRatio'};
     dataOps = {@(x) x, @(x) rad2deg(x), @(x) sqrt(1 - 1./x.^2)};
     dataTitles = {'speed', 'heading change $(\frac{d\theta}{dt})$', 'eccentricity'};
     dataYLabels = {'pixels/s', 'deg/s', ''};
     dataOverlays = {[], {'so.dthetaHiThresh', 'so.dthetaLoThresh'}, {'so.reversalCovThresh', 'so.omegaCovThresh'}};
     overlayLabels = {[], {'d\theta high thresh', 'd\theta low thresh'}, {'reversal min', 'omega max'}};
    
     mirrorOverlays = [false, true, false];
     dataYLimits = {[0 2], [-360, 360], [0 1]};
     dataBackgroundColor = [0 0 0];
     dataAxesColor = [1 1 1];
     dataLineColor = 'w.-';
     dataCurrentMarker = 'w.';
     dataCurrentMarkerSize = 40;
     dataLineWidth = 2;
    
     dataXRange = 30; %time buffer on either side to plot
    
     overlayContour = false;
     overlayMidLine = false;
     interpImage = true;
    
     colormap = gray(256);
    
     delayTime = 0.120; %target interframe interval if playing movie to screen
   
     imAxisSize = 30;
     imBackgroundColor = [0 0 0];
     imAxesColor = [1 1 1];
     imCLim = [10 125];
    
    
     overlayTrack = true;
     locField = 'sloc';
     ptbuffer = 200; %number of points on either side of point to plot on overlay
     trackStyle = 'b-';
     trackWidth = 1;
     ptMarker = 'b.';
     ptMarkerSize = 10;
    
    
     labelRuns = true;
     runColor = 'g';
     labelReorientations = true;
     reorientationColor = 'r';
    
     fontSize = 16;
    
    %worm specific options
     labelSharpTurns = true;
     stbuffer = 0;
    end
    methods(Static)
        function mo = basicAviToDisk(filename) 
            mo = MovieOptions();
            mo.aviName = filename;
            mo.makeAvi = true;
            mo.avioptions = {'FPS', 10, 'COMPRESSION', 'cinepak'};
            mo.stbuffer = 5;
            mo.fontSize = 14;
        end
        function mo = basicOnScreen() 
            mo = MovieOptions();
            mo.imAxesColor = [0 0 0];
            mo.dataAxesColor = mo.imAxesColor;
            mo.dataBackgroundColor = [1 1 1];
            mo.backgroundColor = mo.dataBackgroundColor;
            mo.dataLineColor = 'k-';
            mo.dataCurrentMarker = 'k.';
        end
        function mo = basicImStack(basedir, stub)
            mo = MovieOptions();
            mo.stbuffer = 5;
            mo.fontSize = 14;
            mkdir(basedir);
            mo.jpegstub = fullfile(basedir, stub);
            mo.makeImStack = true;
        end
    end
    
    
end

