classdef MartaExperiment < Experiment
    %An experiment loaded from spine and outline files, instead of .bin
    %files;  for Marta Zlatic's data
    
    properties
    end
    
    methods
        basicCleanup(expt, mintime);
        setSegmentSpeeds (expt, mso);

    end
    
     methods %constructor
        function expt = MartaExperiment(varargin)
            if (~isempty(varargin) && isa(varargin{1}, 'Experiment'))
                fn = intersect(fieldnames(varargin{1}), fieldnames(expt));
                for j = 1:length(fn)
                    expt.(fn{j}) = varargin{1}.(fn{j});
                end
            end
        end
    end
    methods(Static)
        expt = fromFile (dirname);
    end
    
end

