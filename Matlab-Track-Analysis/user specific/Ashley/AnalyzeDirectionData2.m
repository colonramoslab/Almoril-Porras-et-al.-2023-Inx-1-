function [] = AnalyzeDirectionData2 (eset)
%function [] = AnalyzeDirectionData (eset)
%
%This function creates the plots for the Gradient experiments from the expt
%set file.

% for j=1:length(eset)
%     expt(j)=eset(j).expt;
% end

expt=eset.expt;

%initialize variables.
colorcycle = [0 0 1;0 1 0;1 0 0;0.75 0.5 0]';

%first let's store the run variables
tracks=[expt.track];
runs = [tracks.run];

ang = rad2deg([runs.startTheta]);
dur = [runs.runTime];
speed = 66*[runs.pathLength]./[runs.runTime];
angstart=[runs.startTheta];
angend=[runs.endTheta];
dang(1:length(angstart))=0;
for j=1:length(angstart)
    dang(j)=rad2deg(AngleBtwVectors(Ang2Vector(angstart(j)),Ang2Vector(angend(j))));
end

%second let's store the reorientation variables
turns = [tracks.reorientation];
tnHS=[turns.numHS];
pauses=turns(tnHS==0);
turns=turns(tnHS>0);
turn_ang_init = [turns.prevDir];
turn_ang_fin=[turns.nextDir];
turn_dang(1:length(turn_ang_init))=0;
for j=1:length(turn_ang_init)
    turn_dang(j)=AngleBtwVectors(Ang2Vector(turn_ang_init(j)),Ang2Vector(turn_ang_fin(j)));
end

count=1;
for j=1:length(turns)
    if(turns(j).numHS>0)
        max_hs_init(j)=turns(j).headSwing(1).maxTheta;
        max_hs_accept(j)=turns(j).headSwing(end).maxTheta;
    
        newang(j)=AngleAdd(turn_ang_init(j),max_hs_init(j));
    
        if(abs(AngleAdd(newang(j),-pi/2))<abs(AngleAdd(turn_ang_init(j),-pi/2)))
            hs_init_todark(j)=1;
        elseif(abs(AngleAdd(newang(j),-pi/2))>abs(AngleAdd(turn_ang_init(j),-pi/2)))
            hs_init_todark(j)=0;
        else
            hs_init_todark(j)=NaN;
        end
    else
        max_hs_init(j)=NaN;
        max_hs_accept(j)=NaN;
        newang(j)=NaN;
        hs_init_todark(j)=NaN;
    end
    
    for k=1:turns(j).numHS
        turn_hs_ang_init(count)=turn_ang_init(j);
        turn_hs_max(count)=turns(j).headSwing(k).maxTheta;
        newang2(j)=AngleAdd(turn_hs_ang_init(count),turn_hs_max(count));
        if(abs(AngleAdd(newang2(j),-pi/2))<abs(AngleAdd(turn_hs_ang_init(count),-pi/2)))
            turn_hs_todark(count)=1;
        elseif(abs(AngleAdd(newang2(j),-pi/2))>abs(AngleAdd(turn_hs_ang_init(count),-pi/2)))
            turn_hs_todark(count)=0;
        else
            turn_hs_todark(count)=NaN;
        end
        
        if(isnan(turn_hs_todark(count)))
            turn_hs_accept(count)=NaN;
        elseif(k==turns(j).numHS)
            turn_hs_accept(count)=1;
        else
            turn_hs_accept(count)=0;
        end
        count=count+1;
    end
    
%     disp(['Turn init heading=' num2str(rad2deg(turn_ang_init(j)))]);
%     disp(['init head sweep max angle=' num2str(rad2deg(max_hs_init(j)))]);
%     disp(['newang=' num2str(rad2deg(newang(j)))]);
%     disp([ 'oldang add=' num2str(rad2deg(abs(AngleAdd(turn_ang_init(j),pi/2))))]);
%     disp([ 'newang add =' num2str(rad2deg(abs(AngleAdd(newang(j),pi/2))))]);
%     disp(['toward dark=' num2str(hs_init_todark(j))]);
end
turn_ang_init=rad2deg(turn_ang_init);
turn_dang=rad2deg(turn_dang);
max_hs_init=rad2deg(max_hs_init);
max_hs_accept=rad2deg(max_hs_accept);
turn_hs_ang_init=rad2deg(turn_hs_ang_init);
turn_hs_max=rad2deg(turn_hs_max);


%now let's store an array that says which pts are in which quadrants
%and another array that stores them into bins
quad_pts{1}=abs(ang)>=135;
quad_pts{2}=(ang>-135)&(ang<-45);
quad_pts{3}=abs(ang)<=45;
quad_pts{4}=(ang<135)&(ang>45);
quad_turn_pts{1}=abs(turn_ang_init)>=135;
quad_turn_pts{2}=(turn_ang_init>-135)&(turn_ang_init<-45);
quad_turn_pts{3}=abs(turn_ang_init)<=45;
quad_turn_pts{4}=(turn_ang_init<135)&(turn_ang_init>45);
quad_hs_pts{1}=abs(turn_hs_ang_init)>=135;
quad_hs_pts{2}=(turn_hs_ang_init>-135)&(turn_hs_ang_init<-45);
quad_hs_pts{3}=abs(turn_hs_ang_init)<=45;
quad_hs_pts{4}=(turn_hs_ang_init<135)&(turn_hs_ang_init>45);
quad_array=[-180 -90 0 90];
numbins=8;
bin_array = (360/numbins/2-180):(360/numbins):(180-360/numbins/2);
bin_array2 =-180:(360/numbins):180;
for m=1:numbins
    bin_pts{m}=(ang>=bin_array2(m))&(ang<bin_array2(m+1));
end


%okay plot run duration vs. run heading
close all;
figure(1);
hold on;
for q=1:length(quad_pts)
    plot(ang(quad_pts{q}),dur(quad_pts{q}),'.', 'Color',colorcycle(:,q));
end
hold off;
xlabel('Run heading (degrees)')
ylabel('Run duration (s)')
embiggen;
axis([-180 180 0 200]);
figure(20);
hold on;
for m=1:numbins
    v=mod(floor(m/2)+1,5);
    if (v==0); v=1; end
    avg=mean(dur(bin_pts{m}));
    stderr=(std(dur(bin_pts{m}))/sqrt(length(dur(bin_pts{m}))));
    errorbar(bin_array(m),avg,stderr,'.', 'Color',colorcycle(:,v));
end
hold off;
xlabel('Run heading');
ylabel('Run duration (s)');
embiggen;
axis([-190 180 0 40]);



%okay now let's plot run speed vs. run heading
figure(2);
hold on;
for m=1:numbins
    v=mod(floor(m/2)+1,5);
    if (v==0); v=1; end
    avg=mean(speed(bin_pts{m}));
    stderr=(std(speed(bin_pts{m}))/sqrt(length(speed(bin_pts{m}))));
    errorbar(bin_array(m),avg,stderr,'.', 'Color',colorcycle(:,v));
end
hold off;
xlabel('Run heading (degrees)');
ylabel('Run speed (microns/s)');
embiggen;
axis([-190 180 0 300]);

%alright now let's plot the run duration histogram and fit with an
%exponential
figure(3);
dx=5:5:70;
disp(['The mean run time for angles toward -90 deg = ' num2str(mean(dur(quad_pts{2}))) ' s.' ]);
disp(['The mean run time for angles toward 90 deg = ' num2str(mean(dur(quad_pts{4}))) ' s.' ]);
num1=hist(dur(quad_pts{2}),dx);
num2=hist(dur(quad_pts{4}),dx);
% fun=@(x,xdata) x(1)*exp(-xdata./x(2));
% fit1=lsqcurvefit(fun,[50 20], dx, num1);
% fit2=lsqcurvefit(fun,[50 20],dx, num2);
hold on;
plot(dx,num1,'.', 'Color', colorcycle(:,2));
plot(dx,num2,'.', 'Color', colorcycle(:,4));
% plot(dx, fit1(1)*exp(-dx./fit1(2)),'Color', colorcycle(:,1));
% plot(dx,fit2(1)*exp(-dx./fit2(2)),'Color', colorcycle(:,3));
hold off;
drawnow;
xlabel('Run duration (s)');
ylabel('Count');
embiggen;
axis([0 100 0.1 100]);
set(gca,'YScale','log');
% disp(['The fitted tau for angles toward 180 deg = ' num2str(fit1(2)) ' s.' ]);
% disp(['The fitted tau for angles toward 0 deg = ' num2str(fit2(2)) ' s.' ]);




%okay now lets plot the histogram of run heading
figure(4);
num=hist(ang,bin_array2);
tot=sum(num);
num=num/tot;
err=sqrt(num.*(1-num)/tot)*100;
num=num*100;
hold on;
bar(bin_array([1:8]),num([1:8]),'FaceColor', colorcycle(:,1));
bar(bin_array([2 3]),num([2 3]),'FaceColor', colorcycle(:,2));
bar(bin_array([4 5]),num([4 5]),'FaceColor', colorcycle(:,3));
bar(bin_array([6 7]),num([6 7]),'FaceColor', colorcycle(:,4));
errorbar(bin_array,num(1:length(bin_array)),err(1:length(bin_array)),'k.');
hold off;
xlabel('Run heading(deg)');
ylabel('Count (%)');
embiggen;
%axis([-190 180 0 100]);



%plot of run steering
figure(5);
count=1;
quad_lims{1}=[135 225];
quad_lims{2}=[-135 -45];
quad_lims{3}=[-45 45];
quad_lims{4}=[45 135];
% fun2=@(x,xdata) x(1)*exp(-(xdata-x(2)./x(3))^2);
for q=1:length(quad_pts)
    if (count==3) count=5; end
    subplot(4,2,count);
    num=hist(dang(quad_pts{q}), -180:10:180);
    bar(-180:10:180, num/sum(num)*100, 'FaceColor', colorcycle(:,q),'EdgeColor', colorcycle(:,q));
    xlabel('Run heading change(deg)')
    ylabel('Count (%)')
    axis([-180 180 0 20]);
    subplot(4,2,count+2);
    plot(dang(quad_pts{q}),ang(quad_pts{q}),'.', 'Color', colorcycle(:,q));
    hold on;
    if (count==1)
        plot(dang(ang<-135),ang(ang<-135)+360,'.', 'Color', colorcycle(:,q));
    end
    % fit=lsqcurvefit(fun,[50 1 1],dx, num);
    % plot(dx, fit(1)*exp(-(xdata-fit(2)./fit(3))^2),'Color', colorcycle(:,q));
    hold off;
    count=count+1;
    xlabel('Run heading change(deg)');
    ylabel('Run heading (deg)');
    axis([-190 180 quad_lims{q}]);
end

figure(6);
hold on;
for q=1:length(quad_pts)
    errorbar(quad_array(q), mean(DeleteNaNs(dang(quad_pts{q}))),std(DeleteNaNs(dang(quad_pts{q})))/sqrt(length(DeleteNaNs(dang(quad_pts{q})))),'Color', colorcycle(:,q));
end
hold off;
xlabel('Run heading (deg)');
ylabel('Mean Run heading change (deg)');
embiggen;
axis([-190 180 -30 30]);


%plot of angle change in turns
figure(7);
count=1;
for q=1:length(quad_turn_pts)
    if (count==3) count=5; end
    subplot(4,2,count);
    num=hist(turn_dang(quad_turn_pts{q}), -180:10:180);
    bar(-180:10:180, num/sum(num)*100, 'FaceColor', colorcycle(:,q),'EdgeColor', colorcycle(:,q));
    xlabel('Turn heading change(deg)')
    ylabel('Count (%)')
    axis([-190 180 0 20]);
    subplot(4,2,count+2);
    plot(turn_dang(quad_turn_pts{q}),turn_ang_init(quad_turn_pts{q}),'.', 'Color', colorcycle(:,q));
    hold on;
    if (count==1)
        plot(turn_dang(turn_ang_init<-135),turn_ang_init(turn_ang_init<-135)+360,'.', 'Color', colorcycle(:,q));
    end
    % fit=lsqcurvefit(fun,[50 1 1],dx, num);
    % plot(dx, fit(1)*exp(-(xdata-fit(2)./fit(3))^2),'Color', colorcycle(:,q));
    hold off;
    count=count+1;
    xlabel('Turn heading change(deg)');
    ylabel('Turn heading (deg)');
    axis([-190 180 quad_lims{q}]);
end

figure(8);
hold on;
for q=1:length(quad_turn_pts)
    errorbar(quad_array(q), mean(DeleteNaNs(turn_dang(quad_turn_pts{q}))),std(DeleteNaNs(turn_dang(quad_turn_pts{q})))/sqrt(length(DeleteNaNs(turn_dang(quad_turn_pts{q})))),'Color', colorcycle(:,q));
%     plot(turn_ang_init(quad_turn_pts{q}),turn_dang(quad_turn_pts{q}));
%     disp(['mean=' num2str(mean(DeleteNaNs(turn_dang(quad_turn_pts{q})))) ' and sterr=' num2str(std(DeleteNaNs(turn_dang(quad_turn_pts{q})))/sqrt(length(DeleteNaNs(turn_dang(quad_turn_pts{q})))))]);
end
hold off;
xlabel('Turn heading (deg)');
ylabel('Mean Turn heading change (deg)');
embiggen;
axis([-190 180 -30 30]);



%now let's plot the probability of reorientation vs. inst. angle
figure(9);
%first make a histogram of all angles in run
eset=ExperimentSet();
eset.expt=expt;
h=eset.makeHistogram('theta',deg2rad(bin_array2),'run');
re=rad2deg(eset.gatherSubField('run','endTheta'));
bar(bin_array2,hist(re,bin_array2)./h);
xlabel('Turn heading (deg)');
ylabel('Prob. Turn');
embiggen;
% axis([-190 180 0 0.015]);




%plot the histogram of headsweeps per reorientation
figure(10);
bar(0:1:5,hist([turns.numHS], 0:1:5),'k');
xlabel('Num. of head sweeps per turn');
ylabel('Count');
embiggen;
disp(['The mean # of head sweeps = ' num2str(mean([turns.numHS]))]);


%plot the max headsweep angle for each quadrant
figure(11);
count=1;
for q=1:length(quad_turn_pts)
    if (count==3) count=5; end
    subplot(4,2,count);
    hold on;
    num=hist(max_hs_init(quad_turn_pts{q}), -180:10:180);
    bar(-180:10:180, num/sum(num)*100, 'FaceColor', colorcycle(:,q),'EdgeColor', colorcycle(:,q));
    xlabel('Init HS max angle (deg)')
    ylabel('Count (%)')
    axis([-190 180 0 20]);
    subplot(4,2,count+2);
    plot(max_hs_init(quad_turn_pts{q}),turn_ang_init(quad_turn_pts{q}),'.','Color', colorcycle(:,q));
    if (q==1)
        hold on;
        plot(max_hs_init(turn_ang_init<-135), turn_ang_init(turn_ang_init<-135)+360,'.', 'Color', colorcycle(:,q));
        hold off;
    end
    hold off;
    count=count+1;
    ylabel('Turn heading (deg)')
    xlabel('Init. Head Sweep max angle (deg)')
    axis([-190 180 quad_lims{q}]);
end


figure(12);
hold on;
for q=1:length(quad_turn_pts)
    errorbar(quad_array(q), mean(abs(max_hs_init(quad_turn_pts{q}))),std(abs(max_hs_init(quad_turn_pts{q})))/sqrt(length(max_hs_init(quad_turn_pts{q}))),'Color', colorcycle(:,q));
end
hold off;
xlabel('Turn heading (deg)');
ylabel('Init. Head Sweep max angle (deg)');
embiggen;
axis([-190 180 20 80]);

figure(13);
for q=1:length(quad_turn_pts)
    subplot(2,2,q);
    plot(max_hs_accept(quad_turn_pts{q}),turn_ang_init(quad_turn_pts{q}),'.','Color', colorcycle(:,q));
    if (q==1)
        hold on;
        plot(max_hs_accept(turn_ang_init<-135), turn_ang_init(turn_ang_init<-135)+360,'.', 'Color', colorcycle(:,q));
        hold off;
    end
    xlabel('Turn heading (deg)');
    ylabel('Accepted Head Sweep max angle (deg)');
    axis([-190 180 quad_lims{q}]);
end

figure(14);
hold on;
for q=1:length(quad_turn_pts)
    errorbar(quad_array(q), mean(abs(max_hs_accept(quad_turn_pts{q}))),std(abs(max_hs_accept(quad_turn_pts{q})))/sqrt(length(max_hs_accept(quad_turn_pts{q}))),'Color', colorcycle(:,q));
end
hold off;
xlabel('Turn heading (deg)');
ylabel('Accepted Head Sweep max angle (deg)');
embiggen;
axis([-190 180 20 80]);


%what is the probability that the initial headsweep is biased toward dark
figure(15);
perp_pts=quad_turn_pts{1}|quad_turn_pts{3};
%hold on;
% for q=1:length(quad_turn_pts)
%     prob_iniths_todark(q)= sum(hs_init_todark(quad_turn_pts{q}))/length(hs_init_todark(quad_turn_pts{q}));
%     prob_iniths_todark_err(q)=sqrt(prob_iniths_todark(q)*(1-prob_iniths_todark(q))/length(hs_init_todark(quad_turn_pts{q})));
% %     disp([prob_iniths_todark(q)]);
%     errorbar(quad_array(q),prob_iniths_todark(q)*100,prob_iniths_todark_err(q)*100,'.','Color', colorcycle(:,q));
% end
prob_iniths_todark= sum(hs_init_todark(perp_pts))/length(hs_init_todark(perp_pts));
prob_iniths_todark_err=sqrt(prob_iniths_todark*(1-prob_iniths_todark)/length(hs_init_todark(perp_pts)));
disp(['probability init HS toward dark=' num2str(prob_iniths_todark)]);
bar(0:1:1,[prob_iniths_todark*100 (1-prob_iniths_todark)*100],'k');
set(gca, 'XTickLabel', {'Away from Light', 'Toward Light'});
hold on;
errorbar(0:1:1,[prob_iniths_todark*100 (1-prob_iniths_todark)*100],[prob_iniths_todark_err*100 prob_iniths_todark_err*100],'k.');
embiggen;
hold off;
% xlabel('1 = toward Dark,  2 = toward Light')
ylabel('Prob. Init. HS (%)')
% axis([0.5 2.5 0 100]);

%what is the probability that a headsweep toward the dark is accepted
figure(16);
perphs_pts=quad_hs_pts{1}|quad_hs_pts{3};
% hold on;
% for q=1:length(quad_hs_pts)
%     prob_darkhs_accept(q)= sum(DeleteNaNs(turn_hs_todark(quad_hs_pts{q}))&(DeleteNaNs(turn_hs_accept(quad_hs_pts{q}))))/sum(DeleteNaNs(turn_hs_todark(quad_hs_pts{q})));
%     prob_darkhs_accept_err(q)=sqrt(prob_darkhs_accept(q)*(1-prob_darkhs_accept(q))/sum(DeleteNaNs(turn_hs_todark(quad_hs_pts{q}))));
% %     disp([prob_darkhs_accept(q)]);
%     errorbar(quad_array(q),prob_darkhs_accept(q)*100,prob_darkhs_accept_err(q)*100,'.','Color', colorcycle(:,q));
% end
prob_darkhs_accept= sum(DeleteNaNs(turn_hs_todark(perphs_pts))&(DeleteNaNs(turn_hs_accept(perphs_pts))))/sum(DeleteNaNs(turn_hs_todark(perphs_pts)));
prob_darkhs_accept_err=sqrt(prob_darkhs_accept*(1-prob_darkhs_accept)/sum(DeleteNaNs(turn_hs_todark(perphs_pts))));
prob_lighths_accept= sum(DeleteNaNs(~turn_hs_todark(perphs_pts))&(DeleteNaNs(turn_hs_accept(perphs_pts))))/sum(DeleteNaNs(~turn_hs_todark(perphs_pts)));
prob_lighths_accept_err=sqrt(prob_lighths_accept*(1-prob_lighths_accept)/sum(DeleteNaNs(~turn_hs_todark(perphs_pts))));
disp(['prob. away hs accepted=' num2str(prob_darkhs_accept)]);
disp(['prob. toward hs accepted=' num2str(prob_lighths_accept)]);
bar(0:1:1,[prob_darkhs_accept*100 prob_lighths_accept*100],'k');
set(gca, 'XTickLabel', {'Away from Light', 'Toward Light'});
hold on;
errorbar(0:1:1,[prob_darkhs_accept*100 prob_lighths_accept*100],[prob_darkhs_accept_err*100 prob_lighths_accept_err*100],'k.');
embiggen;
hold off;
% xlabel('Turn heading (deg)')
ylabel('Prob. HS accepted (%)')
% axis([-190 180 0 100]);

drawnow;
end