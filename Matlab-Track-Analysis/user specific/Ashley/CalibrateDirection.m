function [dir_mask, cen_pos, dir] = CalibrateDirection()
%
%   
%this program inputs experimental data that is for calibrating the
%directional setup

[FileName, DataPath] = uigetfile('*.*','Select the .bin files you would like to use', 'MultiSelect', 'on');

numfiles=iscellstr(FileName)*length(FileName);
if (numfiles==0) numfiles=1; end

disp('The images will open and line will pop up.');
disp('Drag the line vertices such that it picks the dir of the shadow, press enter.');
disp('Repeat this process for each image.');
for n=1:numfiles
    if (numfiles==1)
        fn=[DataPath FileName];
    else
        fn=[DataPath FileName{n}];
    end
    
    img=imread(fn);
    figure(30);
    imagesc(img);
    colormap('gray');
      
    
    h=impoly(gca,[1000 1000; 1000 1200]);
    setColor(h,[1 0 0]);
    setVerticesDraggable(h,1);
    pause();
    pos=getPosition(h);

    
    shad_init(n,:)=pos(1,:);
    shad_fin(n,:)=pos(2,:);
    cen_pos(n,:)=pos(1,:);
    %773 is 90 um/pix for 2.5 inch post
%     theta(n)=atan2(773,sign(pos(2,2)-pos(1,2))*sqrt((pos(2,1)-pos(1,1))^2+(pos(2,2)-pos(1,2))^2));
    dir(n)=atan2(pos(2,2)-pos(1,2),pos(2,1)-pos(1,1));
    
     
end

num_across=input('How many post positions across are there?');
dir_reshape=mean(reshape(dir,num_across,length(dir)/num_across)');
cen_pos_reshape=mean(reshape(cen_pos(:,1),num_across,length(cen_pos(:,1))/num_across)');
fo=polyfit(cen_pos_reshape,dir_reshape,1);

temp(1:2592)=fo(2)+(fo(1)*[1:1:2592]);
dir_mask=repmat(temp,1944,1);

end

