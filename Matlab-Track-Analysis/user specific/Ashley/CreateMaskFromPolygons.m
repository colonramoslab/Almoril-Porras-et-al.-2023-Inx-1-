function [mask] = CreateMaskFromPolygons (polygon_pos, polygon_value, exclusions, xdimension, ydimension)
%function [mask] = CreateMaskFromPolygons (polygon_pos, polygon_value, exclusions, xdimension, ydimension)
%
%Takes the polygon_pos matrix, which are the position values of the 
%polygons that define the region of interest and makes a mask. The mask is
%not binary, the value is set by the polygon_value, which is an array that
%gives the value for each mask region. exclusions is an array
%that has a 1 to include the polygon and 0 to exclude the polygon from the
%mask. The x and y dimension of the mask is given by xdimension and
%ydimension.

mask=zeros(ydimension,xdimension);
for k=1:round(length(polygon_pos(:,1))/4)
    rectcount=(k-1)*4+1;
    if(exclusions(k)==1)
        mask=MergeMatrix(mask,(poly2mask(polygon_pos(rectcount:(rectcount+3),1),polygon_pos(rectcount:(rectcount+3),2),ydimension,xdimension))*polygon_value(k));
    end
end


% if((blackintop)&&(~IsNumOdd(numacrosstop)))||((~blackintop)&&(IsNumOdd(numacrosstop)))
%         if((k~=ceil(numacrosstop/2))&&(mod(k-ceil(numacrosstop/2),numacrosstop)~=0))
%             rtsidemask=MergeMatrix(rtsidemask,(poly2mask(roipos(rectcount:(rectcount+3),1),roipos(rectcount:(rectcount+3),2),ydimension,xdimension)));
%         end
%     else
%         if(mod(k,numacrosstop)~=0)
%             rtsidemask=MergeMatrix(rtsidemask,(poly2mask(roipos(rectcount:(rectcount+3),1),roipos(rectcount:(rectcount+3),2),ydimension,xdimension)));
%         end
%     end
%     
%     if (blackintop)
%         if(mod(k-ceil(numacrosstop/2),numacrosstop)~=0)
%             ltsidemask=MergeMatrix(ltsidemask,(poly2mask(roipos((rectcount+4):(rectcount+7),1),roipos((rectcount+4):(rectcount+7),2),ydimension,xdimension)));
%         end
%     else
%         if(mod(k,numacrosstop)~=1)
%             ltsidemask=MergeMatrix(ltsidemask,(poly2mask(roipos((rectcount+4):(rectcount+7),1),roipos((rectcount+4):(rectcount+7),2),ydimension,xdimension)));
%         end
%     end
%             
%     if(k>(floor(numacrosstop/2)+(blackintop~=1)))        
%         topmask=MergeMatrix(topmask,(poly2mask(roipos((rectcount+8):(rectcount+11),1),roipos((rectcount+8):(rectcount+11),2),ydimension,xdimension)));
%     end
%     if(k<=(floor(numacrosstop*numacrossside/2)-floor(numacrosstop/2)-(((blackintop~=1)&&(IsNumOdd(numacrossside)))||((blackintop==1)&&(~IsNumOdd(numacrossside))))))
%         botmask=MergeMatrix(botmask,(poly2mask(roipos((rectcount+12):(rectcount+15),1),roipos((rectcount+12):(rectcount+15),2),ydimension,xdimension)));
%     end
%         
% rtsidemask=rtsidemask>0;
% ltsidemask=ltsidemask>0;
% topmask=topmask>0;
% botmask=botmask>0;
% rtsidemask=rtsidemask*1;
% ltsidemask=ltsidemask*2;
% topmask=topmask*3;
% botmask=botmask*4;
% 
% mask=rtsidemask+ltsidemask;
% mask(mask>2)=2;
% mask=mask+topmask;
% mask(mask>3)=3;
% mask=mask+botmask;
% mask(mask>4)=4;
% mask=MergeMatrix(rtsidemask,ltsidemask);
% mask=MergeMatrix(mask,topmask);
% mask=MergeMatrix(mask,botmask);


end