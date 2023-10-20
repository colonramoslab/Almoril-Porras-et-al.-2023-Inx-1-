% %gather all reorientations and runs
% reo=eset.gatherField('reorientation');
% run=eset.gatherField('run');
% 
% %collect number of HS in each reorientation
% reoHS = [reo.numHS];    %creates array of hs numbers where length is total num of reorientation
% 
% %gather light level and derivative of light level for each reorientation
% reodLT = eset.gatherFromSubField('reorientation','dlightTarget','position','mean'); %'position', 'mean' makes this only give one value (the average) for the subfield during all frames of the field
% reoLT = eset.gatherFromSubField('reorientation','lightTarget','position','mean');
% reoLTeti = eset.gatherFromSubField('reorientation','eti','position','mean');
% 
% %gather absolute start frame of all reorientations and runs
% reoStartInd=[reo.startInd]; %Gathers all the start indices for each hs, which is the relative frame number since the beginning of the track
% reoStartFrame = zeros(1,length(reo)); %make matrix of zeros for speed
% for j=1:length(reo)
%     reoStartFrame(j)=reo(j).track.startFrame; %convert the relative startInd to the actual frame number by looking up frame number 
%                                             %stored in each track and adding to startInd
% end
% reoFrameNum=reoStartInd + reoStartFrame;
% reoTimes=eti(reoFrameNum); %takes reoFrameNum for each reorientation and looks up time in eti since they have the same indexing
% 
% %get rid of zeroHS reorientations ie look only at turns
% turnTimes=reoTimes(reoHS>0);   %reoTimesNoZero is an array of all of the times where headsweeps occured
%                                     %the number of elements is equal to the
%                                    %total number of reorientation
%                                     
% %now do same thing for runs
% runStartInd=[run.startInd];
% runStartFrame=zeros(1,length(run));
% for i=1:length(run)
%     runStartFrame(i)=run(i).track.startFrame;
% end
% runFrameNum=runStartInd + runStartFrame;
% runTimes=eti(runFrameNum);
% 
% %gather all light values during runs
% rundLT = eset.gatherField('dlightTarget','run','notlast'); %gather light values during runs except for last run, which is terminated abruptly due to projector shut off
% runLT = eset.gatherField('lightTarget','run','notlast'); %this contains the light values for each frame flagged as a containing a run
% runLTeti=eset.gatherField('eti','run','notlast'); %gather time corresponding to rundLT and runLT
% 
% %gather all target light values on frame by frame basis
% allLT = eset.gatherField('lightTarget');
% alldLT = eset.gatherField('dlightTarget');
% alldLTeti=eset.gatherField('eti');
% 
% %Plot histogram of where turns occur
% fignum = fignum + 1; figure(fignum); clf(fignum);
% hold on
% [y,x]=hist(turnTimes,0:10:endSec);
% plot(x,y);
% title('Hist of turns v. time')
% plot(eti,lt, '-r');
% ylabel('# turns');
% xlabel('time');
% legend('turns','target light');
% axis([0 2500 0 maxLT+2])
% hold off

% %Plot histogram of where pauses occur
% %gather times for zero hs reorientations by looking up when there were no
% %headsweeps
% pauseTimes=reoTimes(reoHS==0);
% %plot pausing hist
% fignum = fignum + 1; figure(fignum); clf(fignum);
% hold on
% [y1,x1]=hist(pauseTimes,0:10:endSec);
% plot(x1,y1);
% title('Hist of pauses v. time')
% plot(eti,lt, '-r');
% ylabel('pauses');
% xlabel('time');
% legend('turning rate','target light');
% axis([0 2500 0 maxLT+2])
% hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    CALCULATE NORMALIZED TURNING  AND PAUSING RATES 
% %    FOR LIGHTS ON/OFF 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %gather all reorientations and runs
% reo=eset.gatherField('reorientation');
% run=eset.gatherField('run');
% 
% %collect number of HS in each reorientation
% reoHS = [reo.numHS];    %creates array of hs numbers where length is total num of reorientation
% 
% %gather light level and derivative of light level for each reorientation
% reodLT = eset.gatherFromSubField('reorientation','dlightTarget','position','mean'); %'position', 'mean' makes this only give one value (the average) for the subfield during all frames of the field
% reoLT = eset.gatherFromSubField('reorientation','lightTarget','position','mean');
% reoLTeti = eset.gatherFromSubField('reorientation','eti','position','mean');
% 
% %collect only turns(throw out reorientations with 0 headsweeps)
% turnsdLT = reodLT(reoHS>0);
% turnsLT = reoLT(reoHS>0);
% turnsLTeti = reoLTeti(reoHS>0);
% 
% %collect only pauses
% pausedLT = reodLT(reoHS==0);
% pauseLT = reoLT(reoHS==0);
% pauseLTeti=reoLTeti(reoHS==0);
% 
% %gather all light values during runs
% rundLT = eset.gatherField('dlightTarget','run','notlast'); %gather light values during runs except for last run, which is terminated abruptly due to projector shut off
% runLT = eset.gatherField('lightTarget','run','notlast'); %this contains the light values for each frame flagged as a containing a run
% runLTeti=eset.gatherField('eti','run','notlast'); %gather time corresponding to rundLT and runLT
% 
% %gather all target light values on frame by frame basis
% allLT = eset.gatherField('lightTarget');
% alldLT = eset.gatherField('dlightTarget');
% alldLTeti=eset.gatherField('eti');
% 
% %calculate conversion factor from frames to minutes
% secPerFrame=eset(1).expt(1).dr.interpTime;
% frame2min=(60/secPerFrame);
% 
% %Gather light values and derivative of light value for pauses
% pauseLT=reoLT(reoHS==0);
% pausedLT=reodLT(reoHS==0);
% 
% 
% %Plot normalized TURNING rate in bins throughout experiment
% bin=1;
% for j=1:(period)/2:round(max(eti))
%     normTurnRate(bin)=(sum((turnsLTeti>j)&(turnsLTeti<j+period/2))/sum((runLTeti>j)&(runLTeti<j+period/2)))*frame2min;
%     bin=bin+1;
% end
%     fignum = fignum + 1; figure(fignum); clf(fignum);
%     plot(normTurnRate);
%     title('Normalized TURNING rate');
%     ylabel('turns per minute');
%     xlabel('bin num');
%     
% %Plot normalized PAUSING rate in bins throughout experiment
% bin=1;
% for j=1:(period)/2:round(max(eti))
%     normPauseRate(bin)=(sum((pauseLTeti>j)&(pauseLTeti<j+period/2))/sum((runLTeti>j)&(runLTeti<j+period/2)))*frame2min;
%     bin=bin+1;
% end
%     fignum = fignum + 1; figure(fignum); clf(fignum);
%     plot(normPauseRate);
%     title('Normalized PAUSING rate');
%     ylabel('pauses per minute');
%     xlabel('bin num');
%     
% %Calculate turning rate index to see if turns go down as a whole throughout
% %expt
% count=1;
% for j=1:2:length(normTurnRate)-1
%     turnIndex(count)=normTurnRate(j)/normTurnRate(j+1);
%     count=count+1;
% end
% 
% fignum = fignum + 1; figure(fignum); clf(fignum);
% plot(turnIndex);
% title('turning rate index');
% ylabel('turns(j)/turns(j+1)');
% xlabel('cycle num of light stimulus');   

%%%%%%%%%%%%%%%%%%% NEW WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %make time range to calc beh metrics
% binLT = lt > 1;
% justOnTimes=zeros(size(binLT));
% for j=1:length(onInds)
%     justOnTimes(onInds(j)+(1:secIntToFrame))=true;
% end
  



%Clac turning rate


% %Plot interval normalized PAUSING rate in bins throughout experiment
% switch(gradientType)
%     case('Square')
%         %sort all vectors for pauses according to time
%         [~,ind2]=sort(pauseLTeti);
%         pauseLTetiS=pauseLTeti(ind2);
%         pauseLTS=pauseLT(ind2);
%         pausedLTS=pausedLT(ind2);
%         
%         %sort all vectors for runs according to time
%         [~,ind2]=sort(runLTeti);
%         runLTetiS=runLTeti(ind2);
%         runLTS=runLT(ind2);
%         rundLTS=rundLT(ind2);
%         
%         %now flag PAUSES that occur at dark->light and light->dark
%         %transitions
%         %changing light value thresh-target light value usually goes from 
%         %0 to 200 in one frame, but sometimes will go in increments of 5ish
%         LTDeltaTh=5;
%         %time span in seconds to acquire pause/run data from after
%         %transition
%         timeSpan=20;
%         
%         %sort PAUSE time and light values into DtoL or LtoD transitions
%         pauseDtoLLT=[];
%         pauseDtoLeti=[];
%         pauseLtoDLT=[];
%         pauseLtoDeti=[];
%         
%         dOnThresh=1;
%         dOffThresh=-1;
%         
%         %sort into dark-->light transitions
%         %will give you 60 frames of transition by using derivative
%         DLInd=find(pausedLTS>dOnThresh);
%         pausedDtoLLT=pausedLTS(DLInd);
%         pauseDtoLeti=pauseLTetiS(DLInd);
% 
%         %find number of pauses that occur when the light is DtoL
%         PauseDtoLNum=sum(pausedLT>dOnThresh);
%         PauseLtoDNum=sum(pausedLT<dOffThresh);
%         
%         
%         
%         LDInd=find(pausedLTS<dOffThresh);
%         
%         
%         for j=timeSpan + 1:length(pauseLTetiS)
%             if pauseLTS(j)-pauseLTS(j-timeSpan)> LTDeltaTh
%                 pauseDtoLLT(end+1)=pauseLTS(j);
%                 pauseDtoLeti(end+1)=pauseLTetiS(j);
%             elseif pauseLTS(j)-pauseLTS(j-timeSpan)< -LTDeltaTh
%                 pauseLtoDLT(end+1)=pauseLTS(j);
%                 pauseLtoDeti(end+1)=pauseLTetiS(j);
%             end
%         end
%         
%         %sort RUN time and light values into DtoL or LtoD transitions
%          runDtoLLT=[];
%          runDtoLeti=[];
%          runLtoDLT=[];
%          runLtoDeti=[];
%          for j=timeSpan + 1:length(runLTetiS)
%             if runLTS(j)-runLTS(j-timeSpan)> LTDeltaTh
%                 runDtoLLT(end+1)=runLTS(j);
%                 runDtoLeti(end+1)=runLTetiS(j);
%             elseif runLTS(j)-runLTS(j-timeSpan)< -LTDeltaTh
%                 runLtoDLT(end+1)=runLTS(j);
%                 runLtoDeti(end+1)=runLTetiS(j);
%             end
%         end      
%         
%         %Get pausing rate for intervals around DARK to LIGHT
%         %transition 
%         bin=1;
%         normPauseIntDtoL=[];
%         normDtoLRun=[];
%         for j=period/2:(period):round(max(eti))
%             normPauseIntDtoL(bin)=sum(pauseDtoLeti>j & pauseDtoLeti<(j+period/2));
%             normDtoLRun(bin)=sum(runDtoLeti>j & runDtoLeti<(j+period/2))*frame2min;
%             bin=bin+1;
%         end
%         
%         bin=1;
%         normPauseIntLtoD=[];
%         normLtoDRun=[];
%         for j=period:period:round(max(eti))
%             normPauseIntLtoD(bin)=sum(pauseLtoDeti>j & pauseLtoDeti<(j+period/2));
%             normLtoDRun(bin)=sum(runLtoDeti>j & runLtoDeti<(j+period/2))*frame2min;
%             bin=bin+1;
%         end
%         
%         %now plot
%         fignum = fignum + 1; figure(fignum); clf(fignum);
%         hold on
%         plot(normPauseIntDtoL,'r.-');   %plot dark to light
%         plot(normPauseIntLtoD,'b.-');   %plot light to dark
%         legend('Dark to light transitions','Light to dark transitions');
%         title('Normalized pausing rates for light transitions')
%         ylabel('pausing rate')
%         xlabel('bin');
%         hold off
%         
%         %calc fraction of pauses that occur at each transitons
%         allPauseTranEti=[pauseDtoLeti pauseLtoDeti];
%         fracPauseDtoL=(length(pauseDtoLeti))/(length(allPauseTranEti));
%         fracPauseLtoD=(length(pauseLtoDeti))/(length(allPauseTranEti));
%         
%         fignum = fignum + 1; figure(fignum); clf(fignum);
%         hold on
%         plot(1,fracPauseDtoL,'r.-','Markersize',22);
%         plot(2,fracPauseLtoD,'b.-','Markersize',22);
%         ylabel('Fraction of pauses at transition');
%         set(gca,'XtickLabel',{'Dark to Light','Light to Dark'});
%         set(gca, 'XTick',[1 2]);
%         title('Normalized pausing rate during light transitions');
%         axis([0 3 0 1]);
%         hold off
%         
% end
%  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate turning rates for WHOLE EXPERIMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%