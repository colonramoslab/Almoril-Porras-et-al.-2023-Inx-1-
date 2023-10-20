classdef MyTest
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        currentState = 'constructed';
    end
    
    methods
        function mt = saveobj(mt)
           disp ('save object called');
           mt.currentState = 'saved';
        end
        
    end
    methods (Static)
        function mt = loadobj(mt)
            disp ('load object called');
            mt.currentState = 'loaded';
        end
    end
    
end

