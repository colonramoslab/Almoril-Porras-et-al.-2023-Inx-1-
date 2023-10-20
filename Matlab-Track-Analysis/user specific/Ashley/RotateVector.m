function [vector] = RotateVector(vector, angle, axis)
%function [vector] = RotateVector(vector, angle, axis)
%
%vector is either a 2D or 3D row vector, angle is the amount
%you want rotated in radians, axis is the axis to rotate about
%(x=1, y=2, z=3, in 2D this isn't used

    %first get size of vector to make sure it is a row vector
    s=size(vector);

    %first make rotation matrix
    rotation=[];
    if(s(1)==1)
        if(s(2)==2)
            rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)];
        elseif(s(2)==3)
            if(axis==1)
                rotation=[1 0 0;0 cos(angle) -sin(angle);0 sin(angle) cos(angle)];
            elseif(axis==2)
                rotation=[cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle)];
            elseif(axis==3)
                rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
            else
                disp('Error in RotateVector, didnt specify axis correctly.');
            end
        else
            disp('Error in RotateVector, vector needs to be either 2D or 3D row vector.');
        end
    else
        disp('Error in RotateVector, vector needs to be either 2D or 3D row vector.');
    end
    
    %then multiply by the rotation matrix
    vector=rotation*vector';
    vector=vector';
        

end