function [] = PlotBoundaryHS(bound_cross)


%initialize variables.


%only look at long boundary crossings and ones where it exits boundary
%properly and not too many headsweeps in one crossing
numframes=[bound_cross.numframes];
exitBD=[bound_cross.exitBD];
enterBD=[bound_cross.enterBD];
nHS=[bound_cross.nHeadsweeps];
bound_cross=bound_cross((numframes>8)&(exitBD==0)&(enterBD==1)&(nHS<6));



%calculate parameters.
% dir=ClosestAng([bound_cross.dir]);
% angin=[bound_cross.angin];
% angin=angin(dir==pi/2);
% bound_dang=[bound_cross.dang];
% bound_dang=bound_dang(dir==pi/2);
% angin=rad2deg(angin);
% bound_dang=rad2deg(bound_dang);
% disp(length(bound_dang))
% clear bound_dang
% clear angin
% clear dir

count=0;
hs_count=0;
for q=1:length(bound_cross)
    if(~isempty(bound_cross(q).decision))
        
        
        
        count=count+1;
        dir(count)=bound_cross(q).dir;
        dir(count)=GetDir(dir(count));
%         disp(length(dir(count)));
%         disp(length(bound_cross(q).angin));
        angin(count)=AngleAdd(bound_cross(q).angin,-1*dir(count));
        angin1(count)=AngleAdd(bound_cross(q).decision(1).heading_stats.init,-1*dir(count));
        bound_dang(count)=bound_cross(q).dang;
        turn_dang(count)=bound_cross(q).decision(1).heading_stats.change;
        hsmax_init(count)=bound_cross(q).decision(1).HS(1).maxTheta;
        hs_accepted_init(count)=bound_cross(q).decision(1).HS(1).accepted;
        if (isnan(hs_accepted_init(count)))
            hs_towardD_init(count)=NaN;
        elseif(sign(angin(count))~=sign(hsmax_init(count)))
            hs_towardD_init(count)=1;
        elseif(sign(angin(count))==sign(hsmax_init(count)))
            hs_towardD_init(count)=0;
        else
            hs_towardD_init(count)=NaN;
            hs_accepted_init(count)=NaN;
        end
        
        
%         disp(['q=' num2str(q) ' and towardD=' num2str(hs_towardD_init(count))]);
%         disp(['angin=' num2str(bound_cross(q).angin)]);
%         disp(['dir=' num2str(bound_cross(q).dir)]);
%         disp(['hsmax=' num2str(hsmax_init(count))]);

        %now do all headsweep props
        for j=1:bound_cross(q).decision(1).nHeadsweeps
            hs_count=hs_count+1;
            hsmax(hs_count)=bound_cross(q).decision(1).HS(j).maxTheta-deg2rad(2);
            hs_accepted(hs_count)=bound_cross(q).decision(1).HS(j).accepted;
            hs_angin(hs_count)=angin(count);
            if(isnan(hs_accepted(hs_count)))
                hs_towardD(hs_count)=NaN;
            elseif(sign(hs_angin(hs_count))~=sign(hsmax(hs_count)))
                hs_towardD(hs_count)=1;
            elseif(sign(hs_angin(hs_count))==sign(hsmax(hs_count)))
                hs_towardD(hs_count)=0;
            else
                hs_towardD(hs_count)=NaN;
                hs_accepted(hs_count)=NaN;
            end
        end
        
        
     
    end
        
        
end


angin=rad2deg(angin);
hsmax_init=rad2deg(hsmax_init);
hsmax=rad2deg(hsmax);
hs_angin=rad2deg(hs_angin);
turn_dang=rad2deg(turn_dang);
bound_dang=rad2deg(bound_dang);
% disp(length(hs_accepted));
% disp(length(DeleteNaNs(hs_accepted)));

%now let's store an array that says which pts are in which quadrants
%and another array that stores them into bins
% quad_pts{1}=abs(angin)>=135;
% quad_pts{2}=(angin>-135)&(angin<-45);
% quad_pts{3}=abs(angin)<=45;
% quad_pts{4}=(angin<135)&(angin>45);
% quad_hs_pts{1}=abs(hs_angin)>=135;
% quad_hs_pts{2}=(hs_angin>-135)&(hs_angin<-45);
% quad_hs_pts{3}=abs(hs_angin)<=45;
% quad_hs_pts{4}=(hs_angin<135)&(hs_angin>45);
% quad_lims{1}=[135 225];
% quad_lims{2}=[-135 -45];
% quad_lims{3}=[-45 45];
% quad_lims{4}=[45 135];
% quad_array=[-180 -90 0 90];
% v=[3 1 4 2 7 5 8 6];
% ang_array=-180:10:180;

quad_pts{1}=abs(angin)>=150;
quad_pts{2}=(angin>-150)&(angin<-90);
quad_pts{3}=(angin>-90)&(angin<-30);
quad_pts{4}=abs(angin)<=30;
quad_pts{5}=(angin<90)&(angin>30);
quad_pts{6}=(angin<150)&(angin>90);
quad_hs_pts{1}=abs(hs_angin)>=150;
quad_hs_pts{2}=(hs_angin>-150)&(hs_angin<-90);
quad_hs_pts{3}=(hs_angin>-90)&(hs_angin<-30);
quad_hs_pts{4}=abs(hs_angin)<=30;
quad_hs_pts{5}=(hs_angin<90)&(hs_angin>30);
quad_hs_pts{6}=(hs_angin<150)&(hs_angin>90);
quad_lims{1}=[150 210];
quad_lims{2}=[-150 -90];
quad_lims{3}=[-90 -30];
quad_lims{4}=[-30 30];
quad_lims{5}=[30 90];
quad_lims{6}=[90 150];
quad_array=[-180 -120 -60 0 60 120];
v=[4 1 5 2 6 3 10 7 11 8 12 9];
ang_array=-180:10:180;
colorcycle = [0 0 1;0 1 0;0.5 0.5 0;1 0 0;0.75 0.5 0.5; 0 0.5 0.5]';


close all;
%plot max hs vs. init ang coming at boundary
figure(1)
for q=1:length(quad_pts)
    subplot(2,3,q);
    plot(hsmax_init(quad_pts{q}),angin(quad_pts{q}),'.','Color', colorcycle(:,q));
    if (q==1)
        hold on;
        plot(hsmax_init(angin<-135), angin(angin<-135)+360,'.', 'Color', colorcycle(:,q));
        hold off;
    end
    xlabel('heading (deg)')
    ylabel('Init. Head Sweep max angle (deg)')
    axis([-180 180 quad_lims{q}]);
end

figure(2)
for q=1:length(quad_hs_pts)
    subplot(2,3,q);
    mean_hsmax(q)=mean(abs(DeleteNaNs(hsmax(quad_hs_pts{q}))));
    err_hsmax(q)=stderr(abs(DeleteNaNs(hsmax(quad_hs_pts{q}))));
    plot(hsmax(quad_hs_pts{q}),hs_angin(quad_hs_pts{q}),'.','Color', colorcycle(:,q));
    if (q==1)
        hold on;
        plot(hsmax(hs_angin<-135), hs_angin(hs_angin<-135)+360,'.', 'Color', colorcycle(:,q));
        hold off;
    end
    xlabel('heading (deg)')
    ylabel('Init. Head Sweep max angle (deg)')
    axis([-180 180 quad_lims{q}]);
end


%what is the probability that the initial headsweep is biased toward dark
figure(3);
hold on;
for q=1:length(quad_pts)
    prob_iniths_todark(q)= sum(hs_towardD_init(quad_pts{q}))/length(DeleteNaNs(hs_towardD_init(quad_pts{q})));
    prob_iniths_todark_err(q)=sqrt(prob_iniths_todark(q)*(1-prob_iniths_todark(q))/length(DeleteNaNs(hs_towardD_init(quad_pts{q}))));
%     disp([prob_iniths_todark(q)]);
    errorbar(quad_array(q),prob_iniths_todark(q)*100,prob_iniths_todark_err(q)*100,'.','Color', colorcycle(:,q));
end
hold off;
xlabel('Turn heading (deg)')
ylabel('Prob. Init. HS to dark (%)')
axis([-180 180 0 100]);

%what is the probability that a headsweep toward the dark is accepted
figure(4);
% disp(~hs_towardD(quad_hs_pts{2}));
% disp(hs_accepted(quad_hs_pts{2}));
hold on;
for q=1:length(quad_hs_pts)
    prob_lighths_accept(q)=sum(DeleteNaNs(~hs_towardD(quad_hs_pts{q}))&(DeleteNaNs(hs_accepted(quad_hs_pts{q}))))/sum(DeleteNaNs(~hs_towardD(quad_hs_pts{q})));
    prob_lighths_accept_err(q)=sqrt(prob_lighths_accept(q)*(1-prob_lighths_accept(q))/sum(DeleteNaNs(~hs_towardD(quad_hs_pts{q}))));
    prob_darkhs_accept(q)= sum(DeleteNaNs(hs_towardD(quad_hs_pts{q}))&(DeleteNaNs(hs_accepted(quad_hs_pts{q}))))/sum(DeleteNaNs(hs_towardD(quad_hs_pts{q})));
    prob_darkhs_accept_err(q)=sqrt(prob_darkhs_accept(q)*(1-prob_darkhs_accept(q))/sum(DeleteNaNs(hs_towardD(quad_hs_pts{q}))));
%     disp([prob_darkhs_accept(q)]);
    errorbar(quad_array(q),prob_darkhs_accept(q)*100,prob_darkhs_accept_err(q)*100,'.','Color', colorcycle(:,q));
    errorbar(quad_array(q),prob_lighths_accept(q)*100,prob_lighths_accept_err(q)*100,'k.');
end
hold off;
xlabel('Turn heading (deg)')
ylabel('Prob. HS to dark accepted (%)')
axis([-180 180 0 100]);


figure(5);
ptstoplot=quad_hs_pts{2}|quad_hs_pts{3};
hold on;
prob_lighths_accept_bar=sum(DeleteNaNs(~hs_towardD(ptstoplot))&(DeleteNaNs(hs_accepted(ptstoplot))))/sum(DeleteNaNs(~hs_towardD(ptstoplot)));
prob_lighths_accept_err_bar=sqrt(prob_lighths_accept(q)*(1-prob_lighths_accept(q))/sum(DeleteNaNs(~hs_towardD(ptstoplot))));
prob_darkhs_accept_bar=sum(DeleteNaNs(hs_towardD(ptstoplot))&(DeleteNaNs(hs_accepted(ptstoplot))))/sum(DeleteNaNs(hs_towardD(ptstoplot)));
prob_darkhs_accept_err_bar=sqrt(prob_darkhs_accept(q)*(1-prob_darkhs_accept(q))/sum(DeleteNaNs(hs_towardD(ptstoplot))));
bar(0:1:1,[prob_darkhs_accept_bar*100 prob_lighths_accept_bar*100],'k');
set(gca, 'XTickLabel', {'','Toward Dark', 'Toward Light',''});
errorbar(0:1:1,[prob_darkhs_accept_bar*100 prob_lighths_accept_bar*100],[prob_darkhs_accept_err_bar*100 prob_lighths_accept_err_bar*100],'k.');  
% disp([prob_darkhs_accept_bar*100 prob_lighths_accept_bar*100]);
hold off;
ylabel('Prob. HS accepted (%)')
embiggen;
axis([-1 2 0 100]);

figure(6);
ptstoplot=quad_pts{2}|quad_pts{3};
hold on;
prob_iniths_todark_bar=sum(DeleteNaNs(hs_towardD_init(ptstoplot)))/length(DeleteNaNs(hs_towardD_init(ptstoplot)));
prob_iniths_todark_err_bar=sqrt(prob_iniths_todark_bar*(1-prob_iniths_todark_bar)/length(DeleteNaNs(hs_towardD_init(ptstoplot))));
prob_iniths_tolight_bar=sum(DeleteNaNs(~hs_towardD_init(ptstoplot)))/length(DeleteNaNs(hs_towardD_init(ptstoplot)));
prob_iniths_tolight_err_bar=sqrt(prob_iniths_tolight_bar*(1-prob_iniths_tolight_bar)/length(DeleteNaNs(hs_towardD_init(ptstoplot))));
bar(0:1:1,[prob_iniths_todark_bar*100 prob_iniths_tolight_bar*100],'k');
set(gca, 'XTickLabel', {'','Toward Dark', 'Toward Light',''});
errorbar(0:1:1,[prob_iniths_todark_bar*100 prob_iniths_tolight_bar*100],[prob_iniths_todark_err_bar*100 prob_iniths_tolight_err_bar*100],'k.');  
% disp([prob_iniths_todark_bar*100 prob_iniths_tolight_bar*100]);
hold off;
ylabel('Prob. Init. HS (%)')
embiggen;
axis([-1 2 0 100]);

figure(7)
errorbar(quad_array, mean_hsmax,err_hsmax,'k.');
xlabel('heading (deg)');
ylabel('Max HS angle (deg)')
embiggen;
axis([-190 180 50 90]);


% figure(8);
% hold on;
% for q=1:length(quad_pts)
%     prob_lighths_accept(q)=sum(DeleteNaNs(~hs_towardD_init(quad_pts{q}))&(DeleteNaNs(hs_accepted_init(quad_pts{q}))))/sum(DeleteNaNs(~hs_towardD_init(quad_pts{q})));
%     prob_lighths_accept_err(q)=sqrt(prob_lighths_accept(q)*(1-prob_lighths_accept(q))/sum(DeleteNaNs(~hs_towardD_init(quad_pts{q}))));
%     prob_darkhs_accept(q)= sum(DeleteNaNs(hs_towardD_init(quad_pts{q}))&(DeleteNaNs(hs_accepted(quad_pts{q}))))/sum(DeleteNaNs(hs_towardD_init(quad_pts{q})));
%     prob_darkhs_accept_err(q)=sqrt(prob_darkhs_accept(q)*(1-prob_darkhs_accept(q))/sum(DeleteNaNs(hs_towardD_init(quad_pts{q}))));
% %     disp([prob_darkhs_accept(q)]);
%     errorbar(quad_array(q),prob_darkhs_accept(q)*100,prob_darkhs_accept_err(q)*100,'.','Color', colorcycle(:,q));
%     errorbar(quad_array(q),prob_lighths_accept(q)*100,prob_lighths_accept_err(q)*100,'k.');
% end
% hold off;
% xlabel('Turn heading (deg)')
% ylabel('Prob. HS to dark accepted (%)')
% axis([-180 180 0 100]);


% %plot turn heading and angle change...
% figure(9)
% count=0;
% for q=1:length(quad_pts)
%     count=count+1;
%     subplot(4,2,v(count));
%     plot(turn_dang(quad_pts{q}),angin(quad_pts{q}),'.','Color', colorcycle(:,q));
%     if (q==1)
%         hold on;
%         plot(turn_dang(angin<-135), angin(angin<-135)+360,'.', 'Color', colorcycle(:,q));
%         hold off;
%     end
%     ylabel('heading (deg)')
%     xlabel('Turn heading change (deg)')
%     axis([-190 180 quad_lims{q}]);
%     
%     count=count+1;
%     subplot(4,2,v(count));
%     bar(ang_array, hist(turn_dang(quad_pts{q}),ang_array),'FaceColor', colorcycle(:,q));
%     xlabel('Turn heading change (deg)')
%     ylabel('Count)')
%     axis([-190 180 0 10]);
% end
% 
% figure(10)
% for q=1:length(quad_pts)
%     mean_turn(q)=mean(turn_dang(quad_pts{q}));
%     err_turn(q)=stderr(turn_dang(quad_pts{q}));
% end
% errorbar([-180 -90 0 90], mean_turn, err_turn,'k.')
% xlabel('heading (deg)')
% ylabel('Turn heading change (deg)')
% axis([-190 180 -30 30]);


%plot boundary heading and angle change...
figure(11);
count=0;
for q=1:length(quad_pts)
    count=count+1;
    subplot(4,3,v(count));
    plot(bound_dang(quad_pts{q}),angin(quad_pts{q}),'.','Color', colorcycle(:,q));
    if (q==1)
        hold on;
        plot(bound_dang(angin<-135), angin(angin<-135)+360,'.', 'Color', colorcycle(:,q));
        hold off;
    end
    ylabel('heading (deg)')
    xlabel('Bound heading change (deg)')
    axis([-190 180 quad_lims{q}]);
    
    count=count+1;
    subplot(4,3,v(count));
    bar(ang_array, hist(bound_dang(quad_pts{q}),ang_array),'FaceColor', colorcycle(:,q));
    xlabel('Bound heading change (deg)')
    ylabel('Count)')
    axis([-190 180 0 10]);
end

figure(12)
for q=1:length(quad_pts)
    mean_bound(q)=mean(bound_dang(quad_pts{q}));
    err_bound(q)=stderr(bound_dang(quad_pts{q}));
end
errorbar(quad_array, mean_bound, err_bound,'k.')
% length(bound_dang)
% mean_bound(1)=mean(bound_dang(abs(angin)>135));
% mean_bound(2)=mean(bound_dang((angin>-135)&(angin<-45)));
% mean_bound(3)=mean(bound_dang(abs(angin)<45));
% mean_bound(4)=mean(bound_dang((angin<135)&(angin>45)));
disp(mean_bound);
xlabel('heading (deg)')
ylabel('Bound heading change (deg)')
axis([-190 180 -60 60]);

end