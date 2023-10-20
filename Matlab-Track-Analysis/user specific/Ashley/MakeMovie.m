function []= MakeMovie(expt,background)
%
%need to comment out hold on drawTrackImage

%initialize avi file
aviobj=avifile('checkbd_pretty','fps',16);

%open figure
close all;
figure(1)


%go through each frame and plot the tracks and maggots
pts=[expt.track.pt];
fileID=fopen(expt.fname);
s_background=size(background);
allblack_background(1:s_background(1),1:s_background(2))=0;

for i=1:4800%length(expt.elapsedTime)
    
    if (mod(i,30)==0) close all; end
%     imagesc(allblack_background)
%     colormap('gray');
%     hold on;
    checkbd_im(1:s_background(1),1:s_background(2))=0;
    
    %figure out pts in frame
    inds=find([pts.ind]==i);
    
    %first plot images
    for j=1:length(inds)
                
        pt_temp=expt.reloadPoint(pts(inds(j)));
        s=size(pt_temp.imData);
        checkbd_im(pts(inds(j)).imOffset(2):(pts(inds(j)).imOffset(2)-1+s(1)),pts(inds(j)).imOffset(1):(pts(inds(j)).imOffset(1)-1+s(2)))=pt_temp.imData;
%         imagesc([pts(inds(j)).imOffset(1) (pts(inds(j)).imOffset(1)+s(2))],[pts(inds(j)).imOffset(2) (pts(inds(j)).imOffset(2)+s(1))], pt_temp.imData)
%         drawnow;
    end
    
%     ih2=imshow(background);
%     colormap('gray');
%     alpha(ih2,0.50);

    checkbd_im=double(checkbd_im)+double(background)*.5;
    imagesc(checkbd_im);
    colormap('gray');
    
    %then plot tracks
    hold on;
    locs=[pts(([pts.ind]>=1)&([pts.ind]<=i)).loc];
    plot(locs(1,:),locs(2,:),'c.','MarkerSize',1);
    hold off;

    
    
    %get frame for avi
    frame=getframe(gcf);
    aviobj = addframe(aviobj,frame);

end

aviobj = close(aviobj);

end
