classdef Simulator < handle
    %Generates a set of simulated tracks in parallel
    %s = Simulator (ntracks, npts) initializes a Simulator with 
    %storage for ntracks & npts
    %
    %set simulation params (depends on type of simulator)
    %
    %Simulator.simulate runs the simulation
    %tracks = Simulator.getTracks outputs a set of tracks from the
    %locations
    
    properties
        ntracks;
        npts;
        loc; %track n, pt p is located at loc(n,:,p);
        timestep = 1; %size of simulation time step
    end
    
    methods %constructor
        function s = Simulator(ntracks, npts)
            s.ntracks = ntracks;
            s.npts = npts;
            s.loc = zeros(ntracks, 2, npts);
        end
    end
    
    methods
        tracks = getTracks(s);
    end
    
    methods (Abstract)
        simulate(s, varargin);
    end
    
end

