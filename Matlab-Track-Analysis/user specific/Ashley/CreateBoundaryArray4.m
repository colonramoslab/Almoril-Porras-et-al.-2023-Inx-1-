function [boundary_im, dir_mask] = CreateBoundaryArray4 (recpos)
%function [boundary_im, dir_mask] = CreateBoundaryArray4 (recpos)
%
%this function creates an array that is the same size as the image. If the
%array has a 0, the point is not at a boundary. If the array has a 1 = right,
%2=left, 3=top, 4=bottom, 5=corner.



%------------------------Initialization---------------------------------
%initialize vars.
cam_y=1944;
cam_x=2592;
width=60; %num of pixels for the width of the region of interest
orientation=[1 0];
size_recpos=size(recpos);
roipos(1:16*floor(size_recpos(1)/4),1:2)=0;
roidir=[];
roiang=[];
roipt=[];

% numacrosstop=input('How many checkers across the top?');
% numacrossside=input('How many checkers across the side?');
% blackintop=input('Is black checker in top left? (Yes=1) (No=0)');
numacrosstop=9;
numacrossside=8;
blackintop=1;




%-------------Make Mask for Polygons that strattle sides-------------------

%now go through the original polygon positions and make new polygons that
%strattle the sides, i.e. the polygon side array, remember in images the 
%y axis is flipped.

for n=1:floor(size_recpos(1)/4)
    topleft=recpos((4*n)-3,:);
    topright=recpos((4*n)-2,:);
    botright=recpos((4*n)-1,:);
    botleft=recpos((4*n),:);
    
    
    %side right rect
    roipos((16*n)-15,:)=[(topright(1)-(width/2)) (topright(2)+(width/2))];
    roipos((16*n)-14,:)=[(topright(1)+(width/2)) (topright(2)+(width/2))];
    roipos((16*n)-13,:)=[(botright(1)+(width/2)) (botright(2)-(width/2))];
    roipos((16*n)-12,:)=[(botright(1)-(width/2)) (botright(2)-(width/2))];
    roidir((4*n)-3,:)=FindBoundaryVector(botright,topright);
    roiang(4*n-3)=AngleBtwVectors(orientation, roidir((4*n)-3,:));
    roipt((4*n)-3,:)=topright;
    
    %side left rect
    roipos((16*n)-11,:)=[(topleft(1)-(width/2)) (topleft(2)+(width/2))];
    roipos((16*n)-10,:)=[(topleft(1)+(width/2)) (topleft(2)+(width/2))];
    roipos((16*n)-9,:)=[(botleft(1)+(width/2)) (botleft(2)-(width/2))];
    roipos((16*n)-8,:)=[(botleft(1)-(width/2)) (botleft(2)-(width/2))];
    roidir((4*n)-2,:)=FindBoundaryVector(topleft,botleft);
    roiang(4*n-2)=AngleBtwVectors(orientation, roidir((4*n)-2,:));
    roipt((4*n)-2,:)=botleft;
    
    %top rect
    roipos((16*n)-7,:)=[(topleft(1)+(width/2)) (topleft(2)+(width/2))];
    roipos((16*n)-6,:)=[(topleft(1)+(width/2)) (topleft(2)-(width/2))];
    roipos((16*n)-5,:)=[(topright(1)-(width/2)) (topright(2)-(width/2))];
    roipos((16*n)-4,:)=[(topright(1)-(width/2)) (topright(2)+(width/2))];
    roidir((4*n)-1,:)=FindBoundaryVector(topright,topleft);
    roiang(4*n-1)=AngleBtwVectors(orientation, roidir((4*n)-1,:));
    roipt((4*n)-1,:)=topleft;
    
    
    %bot rect
    roipos((16*n)-3,:)=[(botleft(1)+(width/2)) (botleft(2)+(width/2))];
    roipos((16*n)-2,:)=[(botleft(1)+(width/2)) (botleft(2)-(width/2))];
    roipos((16*n)-1,:)=[(botright(1)-(width/2)) (botright(2)-(width/2))];
    roipos((16*n),:)=[(botright(1)-(width/2)) (botright(2)+(width/2))];
    roidir((4*n),:)=FindBoundaryVector(botleft,botright);
    roiang(4*n)=AngleBtwVectors(orientation, roidir((4*n),:));
    roipt((4*n),:)=botright;
    
end

for k=1:length(roiang)
    if(roiang(k)==0)
        roiang(k)=0.001;
    end
end


%Calculate the array that excludes polygon sides we don't want.
numpolysides=round(length(roipos(:,1))/4);
numpolys=round(length(roipos(:,1))/16);
exclusions(1:numpolysides)=0;
for k=1:numpolys
    rectcount=(k-1)*4+1; %counts out the number of poly sides
    
    %rt sides
    if((blackintop)&&(~IsNumOdd(numacrosstop)))||((~blackintop)&&(IsNumOdd(numacrosstop)))
        if((k~=ceil(numacrosstop/2))&&(mod(k-ceil(numacrosstop/2),numacrosstop)~=0))
           exclusions(rectcount)=1;
        end
    else
        if(mod(k,numacrosstop)~=0)
            exclusions(rectcount)=1;
        end
    end
    
    %left sides
    if (blackintop)
        if(mod(k-ceil(numacrosstop/2),numacrosstop)~=0)
            exclusions(rectcount+1)=1;
        end
    else
        if(mod(k,numacrosstop)~=1)
            exclusions(rectcount+1)=1;
        end
    end
    
    %tops
    if(k>(floor(numacrosstop/2)+(blackintop~=1)))   
        exclusions(rectcount+2)=1;
    end
    
    %bots
    if(k<=(floor(numacrosstop*numacrossside/2)-floor(numacrosstop/2)-(((blackintop~=1)&&(IsNumOdd(numacrossside)))||((blackintop==1)&&(~IsNumOdd(numacrossside))))))
        exclusions(rectcount+3)=1;
    end
    
end
    
mask=CreateMaskFromPolygons(roipos, exclusions, exclusions, cam_x, cam_y);
dir_mask=CreateMaskFromPolygons(roipos, roiang, exclusions, cam_x, cam_y);




%-----------------------Find corner mask---------------------------------
%calculate the corner positions
size_recpos=size(recpos);
cornpos(1:16*floor(size_recpos(1)/4),1:2)=0;
temp_width=width;
width=width+4;
for n=1:floor(size_recpos(1)/4)
    topleft=recpos((4*n),:);
    topright=recpos((4*n)-1,:);
    botright=recpos((4*n)-2,:);
    botleft=recpos((4*n)-3,:);
    
    
    %top right corner
    cornpos((16*n)-15,:)=[(topright(1)-(width/2)) (topright(2)+(width/2))];
    cornpos((16*n)-14,:)=[(topright(1)-(width/2)) (topright(2)-(width/2))];
    cornpos((16*n)-13,:)=[(topright(1)+(width/2)) (topright(2)-(width/2))];
    cornpos((16*n)-12,:)=[(topright(1)+(width/2)) (topright(2)+(width/2))];
    
    %top left corner
    cornpos((16*n)-11,:)=[(topleft(1)-(width/2)) (topleft(2)+(width/2))];
    cornpos((16*n)-10,:)=[(topleft(1)-(width/2)) (topleft(2)-(width/2))];
    cornpos((16*n)-9,:)=[(topleft(1)+(width/2)) (topleft(2)-(width/2))];
    cornpos((16*n)-8,:)=[(topleft(1)+(width/2)) (topleft(2)+(width/2))];
    
    %bot right corner
    cornpos((16*n)-7,:)=[(botright(1)-(width/2)) (botright(2)+(width/2))];
    cornpos((16*n)-6,:)=[(botright(1)-(width/2)) (botright(2)-(width/2))];
    cornpos((16*n)-5,:)=[(botright(1)+(width/2)) (botright(2)-(width/2))];
    cornpos((16*n)-4,:)=[(botright(1)+(width/2)) (botright(2)+(width/2))];
    
    %bot left corner
    cornpos((16*n)-3,:)=[(botleft(1)-(width/2)) (botleft(2)+(width/2))];
    cornpos((16*n)-2,:)=[(botleft(1)-(width/2)) (botleft(2)-(width/2))];
    cornpos((16*n)-1,:)=[(botleft(1)+(width/2)) (botleft(2)-(width/2))];
    cornpos((16*n),:)=[(botleft(1)+(width/2)) (botleft(2)+(width/2))];

 
end
width=temp_width;

%exclude corners we don't want, i.e. the four of the board
exclusions=[];
exclusions(1:round(length(cornpos(:,1))/4))=0;
for k=1:round(length(cornpos(:,1))/16)
    rectcount=(k-1)*4+1;
    
    %topleft
    if((~blackintop)&&(k==1))
    else
       exclusions(rectcount+3)=1;
    end 
    
    %toprt
    if(k==ceil(numacrosstop/2))
        if((blackintop&&IsNumOdd(numacrosstop))||(~blackintop&&~IsNumOdd(numacrosstop)))
            exclusions(rectcount+2)=1;
        end
    else
        exclusions(rectcount+2)=1;
    end
    
    %botrt
    if(k==ceil(numacrosstop*numacrossside/2))
        if(((~blackintop)&&(IsNumOdd(IsNumOdd(numacrosstop)+IsNumOdd(numacrossside))))||((blackintop)&&(~IsNumOdd(IsNumOdd(numacrosstop)+IsNumOdd(numacrossside)))))
            exclusions(rectcount)=1;
        end
    else
        exclusions(rectcount)=1;
    end
    
    %botleft
    if(k==(ceil(numacrosstop*numacrossside/2)-ceil(numacrosstop/2)+1))
        if((blackintop&&IsNumOdd(numacrossside))||(~blackintop&&~IsNumOdd(numacrossside)))
            exclusions(rectcount+1)=1;
        end
    else
        exclusions(rectcount+1)=1;
    end
end


cornmask=CreateMaskFromPolygons(cornpos, (exclusions*-5), exclusions, cam_x, cam_y);
boundary_im=MergeMatrix(mask, cornmask);

% figure(6)
% imagesc(dir_mask)
% colormap('jet');
% 
% figure(7);
% imagesc(boundary_im);
% colormap('gray');
% hold on;
% for n=1:length(roidir(:,1))
% %     disp(['dir=' num2str(roidir(n,:))]);
% %     disp(['pt=' num2str(roipt(n,:))]);
%     quiver(roipt(n,1),roipt(n,2),roidir(n,1),roidir(n,2),20,'c');
% end
% hold off;
% cd(init.MatlabPath)

TotPix=0;
BndPix=0;
size_bmask=size(boundary_im);
for n=1:size_bmask(1)
    for m=1:size_bmask(2)
        if(boundary_im(n,m)~=0)
            BndPix=BndPix+1;
        end
        TotPix=TotPix+1;
    end
end

% disp(['The % of area covered by boundary = ' num2str(BndPix/TotPix)]);

end