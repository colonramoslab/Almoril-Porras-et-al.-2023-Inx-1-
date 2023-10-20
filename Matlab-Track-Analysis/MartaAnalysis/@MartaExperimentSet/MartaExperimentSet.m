classdef MartaExperimentSet < ExperimentSet
    %extends experiment set to take into account different file loading
    %for marta experiments.
    
    properties
    end
    
    methods (Static)
        eset = fromDirectories(varargin);
        dirlist = findSpineDirectories(varargin);
    end
    
end

