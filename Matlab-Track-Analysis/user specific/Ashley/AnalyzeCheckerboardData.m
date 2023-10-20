function [bound_cross] = AnalyzeCheckerboardData (tracks, checker_im, boundary_im, dir_mask)
%function [bound_cross] = AnalyzeCheckerboardData (tracks, checker_im, boundary_im, dir_mask)
%
%This program first calculates the boundary and then goes through the maggot
%tracks to calculate the parameters inside and outside of the boundary. 
%The parameters for inside the boundary are listed in bound_cross. This
%struct is a list of all of the boundary crossings - started in light/dark/
%other, finished in light/dark/other, num of turns, number of headsweeps,
%incident angle, etc. The parameters for outside of the boundary are listed
%in the spontaneous struct. This struct lists all the turns that occur
%in the interior of the square - number of headsweeps, whether occurs in
%light or dark and all of the runs that occur in the interior of the
%square - duration, velocity, pathlength, etc. and how much time was 
%spent in the interior of the square in the light and dark.


%---------Initialization------------------------------
frame_rate=4;
num_frames_forbndleave=2*frame_rate;


%make sure to add params about headsweep
% MaggotTracks=CalculateInHS(MaggotTracks);

%find boundary
% [boundary_im, dir_mask]=CreateBoundaryArray4(init, polys);
boundary_mask=boundary_im>0;
boundary_mask=boundary_mask+(boundary_im<0);




%go through tracks and if the larva head enters the boundary then mark as
%inBoundary, 1=side, -5=corner
bndcounts=0;
t.inBoundary=[];
t.LDval=[];
for n=1:length(tracks)
%     disp(['track=' num2str(n)]);
    t(n).inBoundary(1:length(tracks(n).pt))=NaN;
    t(n).LDval(1:length(tracks(n).pt))=NaN;
    for p=1:length(tracks(n).pt)
%         disp(['number=' num2str(p)]);
        xpix=round(tracks(n).pt(p).head(1));
        ypix=round(tracks(n).pt(p).head(2));
%         disp(xpix)
%         disp(ypix)
                
        if((ypix<Inf)&&(ypix>500))
            
            if(boundary_mask(ypix,xpix)==1)
                bndcounts=bndcounts+1;
                if(boundary_im(ypix,xpix)==1)
                    t(n).inBoundary(p)=1;
                else
                    t(n).inBoundary(p)=-5;
                end
            else
                t(n).inBoundary(p)=0;
            end
            
            if(checker_im(ypix,xpix)==1)
                t(n).LDval(p)=1;
            elseif(checker_im(ypix,xpix)==0)
                t(n).LDval(p)=0;
            else
                t(n).LDval(p)=NaN;
            end
                
            
        else
            t(n).LDval(p)=NaN;
            t(n).inBoundary(p)=NaN;
        end    
    end  
%      disp(['number=' num2str(p)]);   
end
% disp(['bndcounts=' num2str(bndcounts)]);








%first separate all boundary crossings
bound_cross=[];
Bound_flag=0;
Bound_flag2=0;
numBndCrss=0;
numCount0=0;
ct=0;
time_dark=0;
time_light=0;
for n=1:length(tracks)
%     disp(['track=' num2str(n)]);
    Bound_flag=0;
    Bound_flag2=0;
    numCount0=0;
    for p=1:length(tracks(n).pt)
%         disp(['p=' num2str(p)]);
        if(t(n).inBoundary(p)>0)
                if((Bound_flag==0)&&(Bound_flag2==0))
                    ct=ct+1;
                    Bound_flag=1;
                    Bound_flag2=1;
                    numBndCrss=numBndCrss+1;
                    bound_cross(numBndCrss).track=n;
                    bound_cross(numBndCrss).frame=p;
                    bound_cross(numBndCrss).numframes=1;
                    bound_cross(numBndCrss).enterLD=t(n).LDval(p);
                    bound_cross(numBndCrss).enterBD=t(n).inBoundary(p);
                    bound_cross(numBndCrss).exitLD=-1;
                    bound_cross(numBndCrss).exitBD=-1;
                    numCount0=0;
                elseif((Bound_flag==1)&&(Bound_flag2==1))
                    ct=ct+1;
                    bound_cross(numBndCrss).numframes= bound_cross(numBndCrss).numframes+1;
                    if((p+1)<=length(t(n).LDval))
                        if(t(n).LDval(p+1)<Inf)
                            bound_cross(numBndCrss).exitLD=t(n).LDval(p+1);
                            bound_cross(numBndCrss).exitBD=t(n).inBoundary(p+1);
                        else
                            bound_cross(numBndCrss).exitLD=t(n).LDval(p);
                            bound_cross(numBndCrss).exitBD=-1;
                        end
                            
                    else
                        bound_cross(numBndCrss).exitLD=t(n).LDval(p);
                        bound_cross(numBndCrss).exitBD=-1;
                    end
                    numCount0=0;
                    
                elseif((Bound_flag==1)&&(Bound_flag2==0))
                    
                    
                    ct=ct+1;
                    Bound_flag2=1;
                    bound_cross(numBndCrss).numframes= bound_cross(numBndCrss).numframes+numCount0+1;
                    if((p+1)<=length(t(n).LDval))
                        bound_cross(numBndCrss).exitLD=t(n).LDval(p+1);
                        bound_cross(numBndCrss).exitBD=t(n).inBoundary(p+1);
                    else
                        bound_cross(numBndCrss).exitLD=t(n).LDval(p);
                        bound_cross(numBndCrss).exitBD=-1;
                    end
                    numCount0=0;
                end
        else
                
                numCount0=numCount0+1;
                Bound_flag2=0;
                if((numCount0>=num_frames_forbndleave))
                    Bound_flag=0;
                end
                
        end
            
    end
        
end  
    
% disp(['numBndCrss=' num2str(numBndCrss)]);
numBndCrsspts=0;
for q=1:length(bound_cross)
   numBndCrsspts=numBndCrsspts+bound_cross(q).numframes;
%    disp(bound_cross(q).numframes);
end
% disp(['numBndCrss pts=' num2str(numBndCrsspts)]);
% disp(['ct=' num2str(ct)])



%now go through and figure out Boundary params
%for each Boundary crossing:
%
%1) is it light to dark, dark to light, dark to dark, light to light
%2) angle of incidence entering boundary, angle of incidence leaving
%boundary
%3) # of headsweeps, # of decisions per crossing
%4) for the headsweeps what is the size
% 
%make a new struct... bound_cross.angin = in rad the angle of 
%                     bound_cross.angout =  in rad
%                     bound_cross.enterLD = 1=light, 0=dark for when enter,-5=corner, -1=track ended 
%                     bound_cross.exitLD = same as above except on exit
%                     bound_cross.enterBD = boundary type, -1=track ended 
%                     bound_cross.exitBD = same as above except on exit
%                     bound_cross.time_light = time in seconds
%                     bound_cross.time_dark = time in seconds
%                     bound_cross.nHeadsweeps= num of headsweeps
%                     bound_cross.nDecisions= num of decisions
%                     bound_cross.nRuns= num of runs
%                     bound_cross.headsweep = (decnum, HSnum, max angle, LDval, direction of head at max angle compared to the boundary)
%                     bound_cross.numframes = the number of frames in the
%                       boudary crossing
%                     bound_cross.track = the track the boundary crossing
%                       is in
%                     bound_cross.frame = the pt in the track where the
%                       boundary crossing starts

for q=1:length(bound_cross)
    
    %first find what boundary we are in
    ypix=round(tracks(bound_cross(q).track).pt(bound_cross(q).frame).head(2));
    xpix=round(tracks(bound_cross(q).track).pt(bound_cross(q).frame).head(1));
    bound_cross(q).dir=dir_mask(ypix,xpix);
    
    %next find the angle in and out of the boundary
    frame_et_start=tracks(bound_cross(q).track).pt(bound_cross(q).frame).et;
    frame_et_end=tracks(bound_cross(q).track).pt(bound_cross(q).frame+bound_cross(q).numframes-1).et;
    dq_etis=tracks(bound_cross(q).track).dq.eti;
    [min_diff, dqframe_start]=min(abs(dq_etis-frame_et_start));
    [min_diff, dqframe_end]=min(abs(dq_etis-frame_et_end));
    %do angle in
    if((dqframe_start<=length(dq_etis))&&(dqframe_start>0)) 
        bound_cross(q).angin=tracks(bound_cross(q).track).dq.theta(dqframe_start);
    else
        disp(['dqframe_start=' num2str(dqframe_start) ' and there are only ' num2str(length(dq_etis)) ' dqframes']);
        bound_cross(q).angin=NaN;
    end
    
    %do angle out
    if((dqframe_end<=length(dq_etis))&&(dqframe_end>0)) 
        bound_cross(q).angout=tracks(bound_cross(q).track).dq.theta(dqframe_end);
    else
        disp(['dqframe_end=' num2str(dqframe_end) ' and there are only ' num2str(length(dq_etis)) ' dqframes']);
        bound_cross(q).angout=NaN;
    end 
    
    %plot the change in angle too
    bound_cross(q).dang=AngleAdd(bound_cross(q).angout,-1*bound_cross(q).angin);
    
    %next find time spent in the light and dark
    LD=t(bound_cross(q).track).LDval(bound_cross(q).frame:(bound_cross(q).frame+bound_cross(q).numframes-1));
    bound_cross(q).time_light=length(LD(LD==1))/frame_rate;
    bound_cross(q).time_dark=length(LD(LD==0))/frame_rate;
    
        

    %next find number of decisions and headsweeps per boundary crossing
    Dec_flag=0;
    HS_flag=0;
    matched=0;
    bound_cross(q).nHeadsweeps=0;
    bound_cross(q).nDecisions=0;
%             disp(['track=' num2str(bound_cross(q).track)]);
%             disp(['numframes=' num2str(bound_cross(q).numframes)]);
%             disp(['frame=' num2str(bound_cross(q).frame)]);
     for g=bound_cross(q).frame:(bound_cross(q).frame+bound_cross(q).numframes-1)
            
            %for each pt in the boundary crossing see if in a decision
            matched=0;
            for k=1:length(tracks(bound_cross(q).track).reorientation)
                if((tracks(bound_cross(q).track).reorientation(k).startInd<=g)&&(g<=(tracks(bound_cross(q).track).reorientation(k).endInd)))
                    matched=k;  
                end
            end
            
            %if pt is in a decision and decision hasn't been counted then count the decision
            if ((matched>0)&&(Dec_flag==0))
                bound_cross(q).nDecisions=bound_cross(q).nDecisions+1;
                bound_cross(q).decision(bound_cross(q).nDecisions).heading_stats.init=tracks(bound_cross(q).track).reorientation(matched).prevDir;
                bound_cross(q).decision(bound_cross(q).nDecisions).heading_stats.fin=tracks(bound_cross(q).track).reorientation(matched).nextDir;
                bound_cross(q).decision(bound_cross(q).nDecisions).heading_stats.change=AngleAdd(tracks(bound_cross(q).track).reorientation(matched).nextDir,-1*tracks(bound_cross(q).track).reorientation(matched).prevDir);
%                 bound_cross(q).decision(bound_cross(q).nDecisions).heading_stats.init=ConvertToBoundaryHeading(MaggotTracks(bound_cross(q).track).decision(matched).heading_stats.init, dir_mask, xpix, ypix);
%                 bound_cross(q).decision(bound_cross(q).nDecisions).heading_stats.fin=ConvertToBoundaryHeading(MaggotTracks(bound_cross(q).track).decision(matched).heading_stats.fin, dir_mask, xpix, ypix);
%                 bound_cross(q).decision(bound_cross(q).nDecisions).headin
%                 g_stats.change=ConvertToBoundaryHeading(MaggotTracks(bound_cross(q).track).decision(matched).heading_stats.change, dir_mask, xpix, ypix);
                bound_cross(q).decision(bound_cross(q).nDecisions).HS=tracks(bound_cross(q).track).reorientation(matched).headSwing;
                bound_cross(q).decision(bound_cross(q).nDecisions).nHeadsweeps=tracks(bound_cross(q).track).reorientation(matched).numHS;
                bound_cross(q).nHeadsweeps=bound_cross(q).nHeadsweeps+bound_cross(q).decision(bound_cross(q).nDecisions).nHeadsweeps;
                Dec_flag=matched;
            elseif(Dec_flag==matched)
            else
                Dec_flag=0;
            end
     end
      
end


end
