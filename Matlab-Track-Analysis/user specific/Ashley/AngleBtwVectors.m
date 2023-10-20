function [angle] = AngleBtwVectors(vector1,vector2)
%function [angle] = AngleBtwVectors(vector1,vector2)
%
%This function takes two vectors and finds the angle between them. In 2D, 
%the angle is defined as the counter clockwise rotation to change the direction
%of vector1 into the direction of vector2, so values range from 0 to pi and
%from 0 to -pi. In 3D, the angle is the shortest angle between the two 
%vectors, so values range from 0 to pi. The vectors needn't be normalized.


    if((length(vector1)==2)&&(length(vector2)==2))
        x1=vector1(1);
        y1=vector1(2);
        x2=vector2(1);
        y2=vector2(2);
        angle = mod(atan2(x1*y2-x2*y1,x1*x2+y1*y2),2*pi);
        if(angle>pi)
            angle=angle-2*pi;  
        end
    elseif((length(vector1)==3)&&(length(vector2)==3))
        angle=(atan2(norm(cross(vector2,vector1)),dot(vector2,vector1)));
    else
        disp('Error. Need to input a 2D or 3D vector.');
        angle=NaN;
    end

%     disp(['Angle=' num2str(rad2deg(angle))]);
    
end