classdef BiasedRandomWalkSimulator < Simulator
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        speed = 1; %generic speed
        thetaAxis; %axis over which all reorientation properties are determined
        reoRate; %vs. thetaAxis
        reoMag; %vs. thetaAxis
        biasDirection; %heading direction
        reoBias; %vs. thetaAxis        
        driftCoeff = 0; %diffusion in theta <dtheta^2> = driftCoeff*t;
    end
    
    properties (Access = protected)
        vdir; %velocity angle
    end
    
    methods
        simulate(s, varargin);
        setParams(s, thetaStep, biasDirection, baseRate, rateDelta, baseMag, magDelta, reoBias);
    end
    methods %constructor
        function s = BiasedRandomWalkSimulator(ntracks, npts)
            s = s@Simulator(ntracks,npts);
            s.vdir = 2*pi*rand(1,ntracks);
        end
    end
    
end

