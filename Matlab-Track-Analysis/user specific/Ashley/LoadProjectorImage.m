function [checker_im] = LoadProjectorImage ()
%function [checker_im] = LoadProjectorImage ()
%
%this function loads the projector image and then uses resizable, draggable
%polygons (4-sided) to define the projector image light and dark regions.
%If the projector is all light or all dark there is a separate setting for
%that.


    window_width=0;
    
    %------load the projector image into a matrix--------------------------
    disp('Load the projector image.');
    [FileName,PathName] = uigetfile('*.*','Select the projector image');
    extension=(FileName((length(FileName)-2):length(FileName))); 
    background.image = imread([PathName FileName],extension);
    

    
    
     %-----plot the projector image----------------------------------------
     figure(4);
     colormap(gray);
     imshow(background.image);
     figaxes=get(gcf,'CurrentAxes');
     s=size(background.image);
     
     
     
     
     %----------------Initialize the pattern-----------------------------
    %initialize pattern matrix and fill it depending on the pattern the
    %user selects
    cam_height=s(1);
    cam_width=s(2);
    pattern=zeros(cam_height,cam_width);
    
    
    %find out the lux in the black and white areas
    %set dark to 0 and light to 1
    dklevel=0;
    ltlevel=1;
    background.dklevel=dklevel;
    background.ltlevel=ltlevel;
    pattern=pattern+dklevel;
    
    %find out if checkerboard or all one color
    type=menu('What type of background is it?', 'Checkerboard', 'All one color');
    
    
    
    
    %----------Make Pattern = Checkerboard-------------------------------
    if (type==1) %checkerboard
        
    %either load defined file with polygon positions or create new file.
    background.choice = menu('Do you want to load predefined polygons?','Yes','No'); 
    
    %if load the file then get the 4-sided polygon positions easily
    if background.choice==1  
        disp('Load the struct with the polygon positions.');
        recpos=LoadStruct();
        
    
    %if need to create polygon positions then ask user about the
    %checkerboard and create polygon positions from users answers    
    elseif background.choice==2
        
        %prompt user for # of checkers
        checkerbrd.n=input('How many squares across the top in the checkerboard?');
        checkerbrd.m=input('How many squares down the side in the checkerboard?');
        
        %prompt the user for the corners of a black square on the top edge
        %and a side edge to get the length of the squares and the color of
        %the top left square
        figure(4);
        disp('In Figure 4, select the top corner of the checkerboard');
        disp('on the side with a black square. Then select the bottom');
        disp('corner of the checkerboard on the same side.');
        [board_x,board_y]=ginput(2);
        
        
        %using the users input values go through and make the checkerboard
        %and save the polygon positions that define the checkerboard
        for k=1:length(board_x)
            board_x(k)=round(board_x(k));
            board_y(k)=round(board_y(k));
        end
        background.patternedge_y=board_y(1);
        background.patternedge_x=board_x(1);
        
        checkerbrd.sqlength=round(abs(board_y(2)-board_y(1))/checkerbrd.m);
        checkerbrd.slope=(board_x(2)-board_x(1))/(board_y(2)-board_y(1));
        checkerbrd.theta=-1*atan(checkerbrd.slope);
        if (board_x(1)>(cam_width/2))
            checkerbrd.topleftsq_white=1;
            checkerbrd.topleftsq_x=round((board_x(1)-(checkerbrd.n*checkerbrd.sqlength))*cos(-1*checkerbrd.theta)-(board_y(1)*sin(-1*checkerbrd.theta)));
            checkerbrd.topleftsq_y=round((board_x(1)-(checkerbrd.n*checkerbrd.sqlength))*sin(-1*checkerbrd.theta)+(board_y(1)*cos(-1*checkerbrd.theta))); 
        else
            checkerbrd.topleftsq_white=0;
            checkerbrd.topleftsq_x=board_x(1);
            checkerbrd.topleftsq_y=board_y(1);
        end
        
        count=0;
        rectcount=1;
        for m=1:checkerbrd.m
            for n=1:checkerbrd.n
                count=count+1;
                if(checkerbrd.topleftsq_white&&mod(count,2))
                    recpos(rectcount,:)=[(checkerbrd.topleftsq_x+(checkerbrd.sqlength*(n-1))) (checkerbrd.topleftsq_y+(checkerbrd.sqlength*(m-1)))];
                    recpos(rectcount+1,:)=[(recpos(rectcount,1)+checkerbrd.sqlength) recpos(rectcount,2)];
                    recpos(rectcount+2,:)=[(recpos(rectcount,1)+checkerbrd.sqlength) (recpos(rectcount,2)+checkerbrd.sqlength)];
                    recpos(rectcount+3,:)=[recpos(rectcount,1) (recpos(rectcount,2)+checkerbrd.sqlength)];
                    rectcount=rectcount+4;
                elseif((~checkerbrd.topleftsq_white)&&(~mod(count,2)))
                    recpos(rectcount,:)=[(checkerbrd.topleftsq_x+(checkerbrd.sqlength*(n-1))) (checkerbrd.topleftsq_y+(checkerbrd.sqlength*(m-1)))];
                    recpos(rectcount+1,:)=[(recpos(rectcount,1)+checkerbrd.sqlength) recpos(rectcount,2)];
                    recpos(rectcount+2,:)=[(recpos(rectcount,1)+checkerbrd.sqlength) (recpos(rectcount,2)+checkerbrd.sqlength)];
                    recpos(rectcount+3,:)=[recpos(rectcount,1) (recpos(rectcount,2)+checkerbrd.sqlength)];
                    rectcount=rectcount+4;
                end
                
            end
        end
        
   
        
    else
        disp('ERROR... Did not input choice correctly.');
    end
      
    
    
        
    %Now we can take the polygons that we loaded or made and check their
    %positions if needed
   check=menu('Do you want to check polygon positions?','Yes','No');
   if check==1
       
    %to check just position polygon in figure and set window size to
    %something nice, then ask user to position polygon and press enter in 
    %the matlab window when they are done, finally take position and change
    %array.... 
    for k=1:round(length(recpos(:,1))/4)
        rectcount=(k-1)*4+1;
        figure(4);
        checkerbrd.sqlength=recpos(rectcount+2,1)-recpos(rectcount,1);
        xmin=round(recpos(rectcount,1)-(checkerbrd.sqlength/2)-window_width);
        xmax=round(recpos(rectcount+2,1)+(checkerbrd.sqlength/2)+window_width);
        if xmin<1
            xmin=1;
        end
        if xmax>cam_width
            xmax=cam_width;
        end
        ymin=round(recpos(rectcount,2)-(checkerbrd.sqlength/2)-window_width);
        ymax=round(recpos(rectcount+2,2)+(checkerbrd.sqlength/2)+window_width);
        if ymin<1
            ymin=1;
        end
        if ymax>cam_height
            ymax=cam_height;
        end
        axis([xmin xmax ymin ymax])
        h(k)=impoly(figaxes,recpos(rectcount:(rectcount+3),:));
        setColor(h(k),[1 0 0]);
        doneyet=input('Push enter to set, 0 to exit.');
        if doneyet==0
            break;
        end
        setVerticesDraggable(h(k),0);
        recpos(rectcount:(rectcount+3),:)=getPosition(h(k));
        
    end

    %save new rectangle positions
    SaveStruct(recpos);
   
   end
   

   %now using the polygon positions make the checkerboard mask that
   %perfectly fits the checkerboard.
   mask=zeros(s(1),s(2));
   for k=1:round(length(recpos(:,1))/4)
        rectcount=(k-1)*4+1;
        mask=mask+poly2mask(recpos(rectcount:(rectcount+3),1),recpos(rectcount:(rectcount+3),2),s(1),s(2));
   end
 
   mask=mask>0;
   pattern=(mask*ltlevel);
   pattern=pattern+((pattern<ltlevel)*dklevel);
    
%     %set the background pattern to the lt level if inside a rectangle.
%     for k=1:length(recpos(:,1))
%         for m=round(recpos(k,2)):(round(recpos(k,2))+round(recpos(k,4)))
%             for n=round(recpos(k,1)):(round(recpos(k,1))+round(recpos(k,3)))
%                 pattern(m,n)=ltlevel;
%             end
%         end
%     end




    %----------------Make Pattern = All one Color-------------------------
    %if all one color set everything to either black or white
    elseif (type==2)
        allonecolor=input('Is it all black or all white? (black=0)(white=1)');
        pattern=zeros(cam_height,cam_width);
        if(allonecolor==1)
            pattern=pattern+background.ltlevel;
        elseif(allonecolor==0)
            pattern=pattern+background.dklevel;
        end
    end
    
    
    %-------------Output Pattern-----------------------------------------
    %set the pattern to be output and let user see the pattern in a figure
    background.pattern=pattern;
    
    figure(4)
    imagesc(background.pattern);
    colormap gray;
    colorbar vert;
    
    %plot the new image on top of the old image to see how good it is
%     figure(5);
%     img1=mat2gray(double(background.image));
%     img2=mat2gray(double(background.pattern));
%     ih1=imshow(img1);
%     axis([1 1944 1 2592]);
%     hold on;
%     ih2=imshow(img2);
%     hold off;
%     colormap('gray');
%     alpha(ih2,0.25);
    
    checker_im=background.pattern;
    
    
end
    
    

