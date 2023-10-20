classdef MaggotReorientation < TrackPart
    %groups maggot head swings into reorientations between runs
    
    properties
        numHS = 0;
        headSwing;
        prevDir = NaN;
        nextDir = NaN;
    end
    properties(Dependent = true)
        prevRun;
        nextRun;
    end
    
    
    methods %constructor
        function reo = MaggotReorientation (track, headswings, prevRun, nextRun)
            %if no headswings (headswings = []), then pass prevRun, nextRun
            %to specify interval
            switch nargin
                case 0
                    return;
                case 1
                    reo.track = track;
                    return
                otherwise
                    reo.track = track;
                    reo.headSwing = headswings;
                    existsAndDefault ('prevRun', []);
                    existsAndDefault ('nextRun', []);
                    reo.calculateMetrics (prevRun, nextRun);
            end
        end
    end
    
    methods
        function nr = get.nextRun(obj)
            nr = obj.getAdjacent('next', 'run');
        end
        function pr = get.prevRun(obj)
            pr = obj.getAdjacent('prev', 'run');
        end
    end
    
    
    methods 
        calculateMetrics(reo, prevRun, nextRun);
        draw(reo, varargin)
    end
end

