function [vector] = MakeVector(pt1, pt2)

    %first make a vector from the two points
    vector=pt2-pt1;
    
    
    %now normalize
    vector=vector./norm(vector);
    
    

end