function [data] = AnalyzeTemporalData (eset, dt, ds, dn)
%function [data] = AnalyzeTemporalData (eset, dt, ds, dn)
%
%We want to calculate velocity vs. time and prob. HS vs. time, etc.
%eset is the experiment set. dt is the time between on/off cycles so
%on at 1 s, off at 60 s would be dt=60s. All on/off cycles are condensed to
%one dt*2 second period. ds is the time in seconds that you
%want to mark as stimulated and which starts at time =0 s. dn is the time
%you want to mark as neutral and starts at time = ds. all spontaneous head
%sweeps are labeled as those after time=ds+dn.

    %--------------initialize variables
    frame_rate=4;
    dint=(dt*2*frame_rate); %the length of the interval
    stop_anal=frame_rate*600; %where to stop analyzing, usually at 600 s
    time(1:dint)=(1:dint)/frame_rate;
    avg_vel(1:dint)=0;
    prob_turn(1:dint)=0;
    prob_HS(1:dint)=0;
    num_inturn(1:dint)=0;
    num_inHS(1:dint)=0;
    num_tot(1:dint)=0;
%     close all;

    
      
    
    
    %------Perform Calculations
    %now set the tracks parameters
    tracks=[eset.expt.track];
    turns=eset.gatherField('reorientation');
    HS=eset.gatherSubField('reorientation','headSwing');
    pts=[tracks.pt];
    
    frame=[pts.ind];
    vel=eset.gatherField('speed')*90;
    tim=eset.gatherField('eti');
    turn_ind=[turns.inds];
    HS_ind=[HS.inds];
    
    hs_frame=[HS.startInd];
    hs_max=abs([HS.maxTheta]);
    hs_dur=([HS.endInd]-hs_frame+1)/frame_rate;
    turns_nHS=[turns.numHS];
    count=1;
    hs_iniths=[];
    for j=1:length(turns_nHS)
        temp=[];
        temp(1:turns_nHS(j))=0;
        temp(1)=1;
        %disp(['temp length=' num2str(length(temp))]);
        hs_iniths((end+1):(end+turns_nHS(j)))=temp;
        %disp(hs_iniths)
    end
    hs_endhs=[HS.accepted];
%     disp(length(hs_iniths));
%     disp(length(hs_endhs));
    
%     %take out all data we don't want to analyze
%     frame=frame(frame<=stop_anal);
%     vel=vel(frame<=stop_anal);
%     inturn=inturn(frame<=stop_anal);
%     inHS=inHS(frame<=stop_anal);
%     hs_frame=hs_frame(hs_frame<=stop_anal);
%     hs_max=hs_max(hs_frame<=stop_anal);
%     hs_dur=hs_dur(hs_frame<=stop_anal);
%     hs_iniths=hs_iniths(hs_frame<=stop_anal);
%     hs_endhs=hs_endhs(hs_frame<=stop_anal);
    
    %then reset frame so that everything goes from 1 to dint
    frame=mod(frame,dint);
    frame(frame==0)=dint;
    hs_frame=mod(hs_frame,dint);
    hs_frame(hs_frame==0)=dint;
    turn_ind=mod(turn_ind,dint);
    turn_ind(turn_ind==0)=dint;
    HS_ind=mod(HS_ind,dint);
    HS_ind(HS_ind==0)=dint;
    

    
    %cycle through the frames and set each parameter
    bin_size=40;
    for i=1:dint
        i_bot=i-bin_size;
        i_top=i+bin_size;
        if(i_bot<0) i_bot=1; end
        if(i_top>dint) i_top=dint; end
        
        avg_vel(i)=mean(DeleteNaNs(vel(frame==i)));
        num_inturn(i)=length(turn_ind(turn_ind==i));
        num_inHS(i)=length(HS_ind(HS_ind==i));
        num_tot(i)=sum(frame(frame==i)>0);
        avg_hsmax(i)=mean(hs_max((hs_frame>=i_bot)&(hs_frame<i_top)));
        avg_dur(i)=mean(hs_dur((hs_frame>=i_bot)&(hs_frame<i_top)));
        num_hsend(i)=sum(hs_endhs((hs_frame>=i_bot)&(hs_frame<i_top)))/length(hs_endhs((hs_frame>=i_bot)&(hs_frame<i_top)));
    end
    
    %do some calcs after loop
    temp(1:dint)=0;
    temp(1)=1;
    temp(dint/2)=-1;
    LDval=cumsum(temp);
%     LDval=0;
           
    prob_turn=num_inturn./num_tot;
    prob_HS=num_inHS./num_tot;
    
    
        %Now we also want to find spontaneous and stimulated headsweeps and
    %see which is bigger
    %let's mark frames as spontaneous or stimulated
    dstim= frame_rate*ds;  % first 2 seconds are stim
    dneutral = frame_rate*dn; %next 5 seconds are neutral 
    stim_start=0:dint:(stop_anal-dint);
    stim_start(1)=1;
        

    %set an array that is 1 when stimulated
    temp(1:4800)=0;
    temp(stim_start) = 1;
    temp(stim_start+dstim) = -1;
    stim = cumsum(temp);
    %now make an array that is 1 when neutral
    temp(1:4800)=0;
    temp(stim_start+dstim) = 1;
    temp(stim_start+dstim+dneutral) = -1;
    neutral=cumsum(temp);
    %finally make an array that is 1 when spontaneous
    temp(1:4800)=0;
    temp(stim_start+dstim+dneutral) = 1;
    temp(stim_start+dint) = -1;
    spon=cumsum(temp);
    
    %now save all the stimulated hs max angles and spont. max angles
    %together... head sweep has to start in the region to be counted.
    stimhs_max=hs_max((stim(hs_frame)==1)&(hs_iniths==1));
    sponhs_max=hs_max((spon(hs_frame)==1)&(hs_iniths==1));
    neuths_max=hs_max((neutral(hs_frame)==1)&(hs_iniths==1));
    
    stimhs_dur=hs_dur((stim(hs_frame)==1)&(hs_iniths==1));
    sponhs_dur=hs_dur((spon(hs_frame)==1)&(hs_iniths==1));
    neuths_dur=hs_dur((neutral(hs_frame)==1)&(hs_iniths==1));
    
    stimhs_endhs=hs_endhs((stim(hs_frame)==1)&(hs_iniths==1));
    sponhs_endhs=hs_endhs((spon(hs_frame)==1)&(hs_iniths==1));
    neuths_endhs=hs_endhs((neutral(hs_frame)==1)&(hs_iniths==1));
  
    
    

    %-------------store data
    data.time=time;
    data.LDval=LDval;
    data.vel=avg_vel;
    data.prob_turn=prob_turn;
    data.prob_HS=prob_HS;
    data.avg_hsmax=rad2deg(avg_hsmax);
    data.avg_dur=avg_dur;
    data.prob_hsend=num_hsend;
    
    
    %--------------plot the figures   
    figure(1);
    subplot(4,1,1);
    plot(data.time,data.LDval, 'b');
    xlabel('Time (s)');
    ylabel('Light on?');
    axis([0 dt*2 -0.1 1.1]);
        
    subplot(4,1,2);
    vel=eset.gatherField('speed');
    tim=eset.gatherField('eti');
    vel_bin=1;
    [x,meany,standarderror,standarddeviation] = meanyvsx (tim, vel, 0:vel_bin:1200);
    for j=1:10
        start_j=1+((j-1)*dt*2/vel_bin);
        end_j=start_j+(dt*2/vel_bin)-1;
        meany_avg(j,:)=meany(start_j:end_j);
    end
    x=vel_bin:vel_bin:(dt*2);
    meany=mean(meany_avg);
    plot(x,meany, 'b');
    xlabel('Time (s)');
    ylabel('Velocity (microns/s)');
    axis([0 dt*2 4 6]);
    
    subplot(4,1,3);
    plot(data.time,data.prob_turn, 'b');
    xlabel('Time (s)');
    ylabel('Prob. turn');
    axis([0 dt*2 0 .6]);
    
    subplot(4,1,4);
    plot(data.time,data.prob_HS, 'b');
    xlabel('Time (s)');
    ylabel('Prob. HS');
    axis([0 dt*2 0 .6]);
    
    
% %     make the vertical lines
% %     vert_lines_y=repmat(stim_start,1,2);
% %     vert_lines_y=sort(vert_lines_y);
% %     vert_lines_x=repmat([0 1],1,length(stim_start));
% %     vert_lines_x2=repmat([50 150],1,length(stim_start));
%     for j=1:length(stim_start)
%         subplot(4,1,1);
%         hold on;
%         plot([stim_start(j)/frame_rate stim_start(j)/frame_rate],[0 1],'r');
%         hold off;
%         subplot(4,1,2);
%         hold on;
%         plot([stim_start(j)/frame_rate stim_start(j)/frame_rate],[50 150],'r');
%         hold off;
%         subplot(4,1,3);
%         hold on;
%         plot([stim_start(j)/frame_rate stim_start(j)/frame_rate],[0 1],'r');
%         hold off;
%         subplot(4,1,4);
%         hold on;
%         plot([stim_start(j)/frame_rate stim_start(j)/frame_rate],[0 1],'r');
%         hold off;
%     end

    
    %plot the headsweep variables
    figure(2);
    subplot(4,1,1);
    plot(data.time,data.LDval, 'b');
    xlabel('Time (s)');
    ylabel('Light on?');
    axis([0 dt*2 -0.1 1.1]);
    
    subplot(4,1,2);
    plot(data.time,data.avg_hsmax, 'b');
    xlabel('Time (s)');
    ylabel('Max HS ang (deg)');
    axis([0 dt*2 60 80]);
    
    subplot(4,1,3);
    plot(data.time,data.avg_dur, 'b');
    xlabel('Time (s)');
    ylabel('Duration (s)');
    axis([0 dt*2 1 2]);
    
    subplot(4,1,4);
    plot(data.time,data.prob_hsend, 'b');
    xlabel('Time (s)');
    ylabel('Prob. accept headsweep');
    axis([0 dt*2 0.5 0.8]);


    

    
    
    
    %head sweep angle
    figure(3);
    hsmax_avg=rad2deg([mean(abs(stimhs_max)) mean(abs(sponhs_max)) mean(abs(neuths_max))]);
    hsmax_err=rad2deg([stderr(abs(stimhs_max)) stderr(abs(sponhs_max)) stderr(abs(neuths_max))]);
    bar([0 1 2],hsmax_avg,'k');
    set(gca, 'XTickLabel', {'Stim', 'Spon', 'Neutral'});
    hold on;
    errorbar([0 1 2],hsmax_avg,hsmax_err,'k.');
    embiggen;
    hold off;
    ylabel('mean HS max angle');
    
    %head sweep duration
    figure(4)
    hsdur_avg=[mean(stimhs_dur) mean(sponhs_dur) mean(neuths_dur)];
    hsdur_err=[stderr(stimhs_dur) stderr(sponhs_dur) stderr(neuths_dur)];
    bar([0 1 2],hsdur_avg,'k');
    set(gca, 'XTickLabel', {'Stim', 'Spon', 'Neutral'});
    hold on;
    errorbar([0 1 2],hsdur_avg,hsdur_err,'k.');
    embiggen;
    hold off;
    ylabel('mean HS duration');
    
    %head sweep rejection
    figure(5)
    hsend_num=[length(stimhs_endhs) length(sponhs_endhs) length(neuths_endhs)];
    hsend_avg=[sum(stimhs_endhs) sum(sponhs_endhs) sum(neuths_endhs)]./hsend_num;
    hsend_err=sqrt(hsend_avg.*(1-hsend_avg)./hsend_num);
    bar([0 1 2],hsend_avg,'k');
    set(gca, 'XTickLabel', {'Stim', 'Spon', 'Neutral'});
    hold on;
    errorbar([0 1 2],hsend_avg,hsend_err,'k.');
    embiggen;
    hold off;
    ylabel('Prob. HS acceptance');

end