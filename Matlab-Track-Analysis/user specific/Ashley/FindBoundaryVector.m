function [vector] = FindBoundaryVector(pt1, pt2)
%function [vector] = FindBoundaryVector(pt1, pt2)
%
%pt1 and pt2 are 2D cartesian coordinates (row vectors)
%input the points such that pt2-pt1 creates a vector that points in a
%direction that is a clockwise 90 degree rotation (-90) from pointing to 
%the dark square:  
%
%dark
%-->
%light
   
    %make a vector from the points and rotate it 90 degrees
    vector=MakeVector(pt1,pt2);
    vector=RotateVector(vector,(pi/2), 3);
    
    
    
    


end