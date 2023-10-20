addpath('basic routines')

display(' ________________________________')
display('|                                |')
display('| WELCOME TO TRACK ANALYSIS 2.3! |')
display('|   (developed by Daniel Kim!)   |')
display('|________________________________|')

disp(' ')
disp(' ')
disp('Let''s load your BIN files!')

disp(' ')
exs = ExperimentSet.fromFiles('minpts',200);

maxframe = 3;
maxdist = 4.5;
minspeed = 0.4;
mindist = 80;
minpts = 400;
prunemarg = 200;
trimpad = 0;

line='~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
line2='\n-------------------------------\n\n';

clean = ESetCleaner();

inds = findstr(exs.expt(1).fname,'\');
ind = inds(numel(inds));

disp(' ')
disp('ALL EXPERIMENTS LOADED.')
disp('LET''S GET THIS PARTY STARTED.')

multi=false; 
for i=1:numel(exs.expt)+1
    clf;
    close all;
    satisfied = false;
    segmented = false;
    disp(' ')
    disp(line)
    disp(' ')
    tempset = ExperimentSet();
    if multi
        reply = input('Would you like to perform a group analysis?\n  0. No\n  1. Yes\nYour Response: ');
        if reply == 0
            satisfied=true;
        else
            tempset = exs;
        end
    else
        ex = exs.expt(i);
        tempset.expt=[ex];
        
        
        flen = numel(ex.fname);
        open([ex.fname(1:flen-4) '_header.txt']);
        
        disp(['Experiment ' num2str(i) ': ' ex.fname(ind+1:numel(ex.fname))])
        disp(' ')
        disp('  Analysis Rectangle')
        x0 = input('  | x0 = ');
        y0 = input('  | y0 = ');
        x1 = input('  | x1 = ');
        y1 = input('  | y1 = ');
        
        disp(' ')
        disp('  Foreground Image')
        fgim = input('  | dir = ','s');
        
        t0 = ex.elapsedTime(1);
        t1 = ex.elapsedTime(numel(ex.elapsedTime));
    end
    
    if ~satisfied
        disp(' ')
        disp('Locate the directory to save figures:')
        figsdir = uigetdir('Results');
        disp([' ' figsdir]);
    end
    
    while ~satisfied
        figs2save=[];
        fignames={};
        reply = input([line2 'What would you like to now?\n  0. Next Experiment\n  1. View Tracks\n  2. Stitch Tracks\n  3. Clean Tracks\n  4. Prune Tracks\n  5. Trim Tracks\n  6. Segment Tracks\n  7. Plot CoM\n  8. Experiment Info\n  9. Clear/Close Figures\n  *. Reload Experiment\nYour Response: '],'s');
        if reply == '*'
            num = -1;
        else
            num = str2num(reply);
        end
        if num == 0
            disp(' ')
            disp('NEXT')
            satisfied=true;
            if i==numel(exs.expt)
                multi=true;
            end
        
        elseif num == 1
            disp(' ')
            disp('VIEW ''EM')
            figs2save(length(figs2save)+1)=figure;
            fignames(length(figs2save))={'tracks_all'};
            alltracks = [tempset.expt.track];
            alltracks.plotPath('displacement');
            hold on
            plot (0,0, 'r.', 'MarkerSize', 20);
            hold off
            
            if ~multi
                if ~isempty(fgim)
                    figure;
                    image(imread(fgim,'JPEG'))
                end
                disp(' ')
                reply = input('  View Individual Tracks?\n  | 0. No\n  | 1. Yes\n  | Your Response: ');
                if reply==1
                    satisfied2 = false;
                    while ~satisfied2
                        disp(' ')
                        disp('  Track Listings:');
                        disp('  | 0. I''m done. Take me home.');
                        spds = ex.gatherField('speed','mean');
                        for tnum=1:length(ex.track)
                            tr=ex.track(tnum);
                            disp(['  | ' num2str(tnum) '. Track ' num2str(tnum) ': ' num2str(spds(tnum)) ' px/s  ' num2str(tr.npts) ' points'])
                        end
                         disp('  | *. All of ''em.');
                        trs=input('  | Desired Tracks: ','s');
                        if trs == '0'
                            satisfied2 = true;
                        elseif trs == '*'
                            figs2save(length(figs2save)+1)=figure;
                            fignames(length(figs2save))={'tracks'};
                            ex.track(1).plotPath;
                            hold on
                            for q=1:length(ex.track)-1
                                ex.track(q+1).plotPath;
                            end
                            hold off;
                        else
                            inds = findstr(trs,',');
                            inds(length(inds)+1)=length(trs)+1;
                            figs2save(length(figs2save)+1)=figure;
                            fignames(length(figs2save))={'tracks'};
                            tind=str2num(trs(1:inds(1)-1)); %#ok<*ST2NM>
                            ex.track(tind).plotPath;
                            disp(' ')
                            disp(['    [ Track ' num2str(tind) ': ' num2str(spds(tind)) ' px/s  ' num2str(ex.track(tind).npts) ' points'])
                            hold on
                            for q=1:length(inds)-1
                                tind=str2num(trs(inds(q)+1:inds(q+1)-1)); %#ok<*ST2NM>
                                ex.track(tind).plotPath;
                                disp(['    [ Track ' num2str(tind) ': ' num2str(spds(tind)) ' px/s  ' num2str(ex.track(tind).npts) ' points'])
                            end
                            
                            hold off
                        end
                    end
                end
            end

        elseif num == 2
            disp(' ')
            reply = input('  Stitch Params (1: new, ENTER: default/last) ');
            if reply==1
                maxframe=input('  |    maxFrame = ');
                maxdist=input('  | maxDistance = ');
            else
                disp(['  |    maxFrame = ' num2str(maxframe)]);
                disp(['  | maxDistance = ' num2str(maxdist)]);
            end
            disp(' ')
            disp('STITCH ''EM')
            disp(' ')
            tempset.executeExperimentFunction('stitchTracks',maxframe,maxdist);
        
        elseif num == 3
            disp(' ')
            reply = input('  Clean Params (1: new, ENTER: default/last) ');
            if reply==1
                minspeed=input('  | minSpeed = ');
                mindist=input('  |  minDist = ');
                minpts=input('  |   minPts = ');
            else
                disp(['  | minSpeed = ' num2str(minspeed)]);
                disp(['  |  minDist = ' num2str(mindist)]);
                disp(['  |   minPts = ' num2str(minpts)]);
            end
            clean.minSpeed = minspeed;
            clean.minDist = mindist;
            clean.minPts = minpts;
            
            disp(' ')
            disp('REPORT:')
            clean.getReport(tempset);
            
            disp(' ')
            reply=input('Keep cleaning tracks (y/n)? ','s');
            if reply=='y'
                clf;
                close all;
                disp('CLEAN ''EM')
                disp(' ')
                clean.clean(tempset);
            else
                disp('DON''T CLEAN ''EM')
            end 
                
        elseif num == 4
            if ~multi
                disp(' ')
                disp('PRUNE ''EM')
                disp(' ')
                tracknum = numel(ex.track);
                reply = input('  Prune Params (1: new, ENTER: default/last) ');
                if reply==1
                    prunemarg=input('  | pruneMargin = ');
                else
                    disp(['  | pruneMargin = ' num2str(prunemarg)]);
                end
                
                for q=1:length(ex.track)
                    starts(q)=ex.track(q).startFrame();
                end
                starts=sort(starts);
                gotit=false;
                q=1;
                while ~gotit && q<=length(starts)
                    yeah=0;
                    for q2=1:length(starts)
                        if starts(q)>=ex.track(q2).startFrame() && starts(q)<=ex.track(q2).endFrame()
                            yeah=yeah+1;
                        end
                    end
                    yeahs(q)=yeah/length(ex.track);
                    if yeah/length(ex.track) > 0.33
                        gotit=true;
                    else
                        q=q+1;
                    end
                end
                
                if ~gotit
                    [m,q]=max(yeahs); %#ok<*ASGLU>
                end
                timeAxis = t0:1:t1;
                [time, centerofmass] = tempset.meanField2vsField1('eti', 'sloc', timeAxis);
                center = centerofmass(1,starts(q)/2);
                
                ex.pruneTracks([t0 t1],[center-prunemarg y0 center+prunemarg y1]);
                disp(' ')
                disp([num2str(tracknum - numel(ex.track)) ' tracks pruned'])
            else
                disp(' ')
                disp('Pruning not allowed for group analysis!')
                disp(' ')
            end
        
        elseif num == 5
            if ~multi
                disp(' ')
                disp('TRIM ''EM')
                disp(' ')
                tracknum = numel(ex.track);
                reply = input('  Trim Params (1: new, ENTER: default/last) ');
                if reply==1
                    trimpad=input('  | trimPadding = ');
                    t2=input('  |          t0 = ');
                    if t2>t0
                        t0=t2;
                    end
                    t2=input('  |          t1 = ');
                    if t2<t1
                        t1=t2;
                    end
                else
                    disp(['  | trimPadding = ' num2str(trimpad)]);
                    disp(['  |          t0 = ' num2str(t0)]);
                    disp(['  |          t1 = ' num2str(t1)]);
                end
                ex.trimTracks([t0 t1],[x0+trimpad y0+trimpad x1-trimpad y1-trimpad]);
                disp(' ')
                disp([num2str(tracknum - numel(ex.track)) ' tracks removed'])
            else
                disp(' ')
                disp('Trimming not allowed for group analysis!')
                disp(' ')
            end
        
        elseif num == 6
            disp(' ')
            disp('SEGMENT ''EM')
            wso = WormSegmentOptions;
            if ~segmented
                tempset.executeExperimentFunction('segmentTracks', wso);
                segmented = true;
                disp('-> FINITO.');
            else
               disp('-> ALREADY FINITO.'); 
            end
            satisfied2 = false;
            while ~satisfied2
                disp(' ')
                reply = input('  Further Analysis?\n  | 0. My work here is done.\n  | 1. HISTOGRAM: heading angle v. time (all)\n  | 2. HISTOGRAM: heading angle v. time (runs)\n  | 3. HISTOGRAM: reorientation rate v. heading angle\n  | 4. SCATTER: angle change v. starting angle\n  | 5. AUTO-CORRELATE: path direction\n  | 6. AUTO-CORRELATE: change in heading over time\n  | 7. CROSS-CORRELATE\n  | 8. TEST: random noise\n  | 9. Clear/Close Figures\n  | Your Response: ');
                tempset.defaultTitle = 'Segmentation Analysis';
                thetaAxis = deg2rad(0:15:345);
                if reply==1
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'theta_all'};
                    tempset.makeHistogram('theta', thetaAxis, 'polar', true, 'r2d', true);
                    snapnow;
                elseif reply==2
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'theta'};
                    tempset.makeHistogram('theta', thetaAxis, 'runs','polar', true, 'r2d', true);
                    snapnow;
                elseif reply==3
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'reorientation'};
                    tempset.makeReorientationHistogram('theta', thetaAxis, 'polar', true, 'r2d', true);
                    snapnow;
                elseif reply==4
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'dtheta_scatter'};
                    runstart = tempset.gatherSubField('run', 'startTheta');
                    runend = tempset.gatherSubField('run', 'endTheta');
                    dt = diff(unwrap([runstart;runend]));
                    runstart = mod(rad2deg(runstart), 360);
                    dt = rad2deg(dt);
                    [rs, meanchange, stderrchange] = meanyvsx(runstart, dt, 0:30:360);
                    plot (runstart, dt, 'k.'); hold on
                    errorbar(rs, meanchange, stderrchange, 'r', 'LineWidth', 2); hold off
                    snapnow;
                elseif reply==5
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'ac_direction'};
                    [ac, np, tx] = tempset.autocorrelate('vnorm');
                    [acr, npr, txr] = tempset.autocorrelate('vnorm','withinRuns',true);
                    semilogy(tx(ac>0), ac(ac>0)./np(ac>0), 'b.',txr(acr>0), acr(acr>0)./npr(acr>0),'g.');
                    xlim([0 600]);
                    xlabel ('$\tau$ (s)','Interpreter', 'Latex');
                    ylabel('$\langle\hat{v}(t)\cdot\hat{v}(t + \tau)\rangle$','Interpreter','Latex');
                    title ('Auto-Correlation of velocity direction');
                    legend('over whole track', 'within runs');
                    ylim([0.01 1])
                    snapnow
                elseif reply==6
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'ac_heading'}; %#ok<*SAGROW>
                    [ac, np, tx] = tempset.autocorrelate('deltatheta','withinRuns',true);
                    plot(tx, ac./np);
                    xlim([0 100])
                    xlabel ('$\tau$ (s)','Interpreter', 'Latex');
                    ylabel('$\langle\dot{\theta}(t)\ast\dot{\theta}(t + \tau)\rangle$','Interpreter','Latex');
                    title ('Auto-Correlation of heading angle change');
                    snapnow
                elseif reply==7
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'cc'};
                    [xc, np, tx] = tempset.crosscorrelate('deltatheta', 'ddtheta', 'withinRuns', true);
                    plot (tx, xc./np); xlim([-50 50])
                    xlabel ('$\tau$ (s)','Interpreter', 'Latex');
                    ylabel('$\langle\ddot{\theta}(t)\ast\dot{\theta}(t + \tau)\rangle$','Interpreter','Latex');
                    title ('Cross correlation of first and second derivatives of heading');
                    snapnow
                elseif reply==8
                    figs2save(length(figs2save)+1)=figure;
                    fignames(length(figs2save))={'error'};
                    t = tempset.expt(1).track(1);
                    t.dq.randomcrap = randn(size(t.dq.eti));
                    t.dq.srandomcrap = lowpass1D(t.dq.randomcrap, t.dr.smoothTime/t.dr.interpTime);
                    t.dq.dsrandomcrap = deriv(t.dq.srandomcrap, t.dr.derivTime/t.dr.interpTime);
                    [xccrap,npcrap,txcrap] = t.crosscorrelate('srandomcrap', 'dsrandomcrap');
                    
                    dgc = -conv(gausskernel(t.dr.smoothTime/t.dr.interpTime), dgausskernel(t.dr.derivTime/t.dr.interpTime),'same');
                    dgc = dgc*max(xccrap(round(length(xccrap)/2) + [-100:100])./npcrap(round(length(xccrap)/2) + [-100:100]))/max(dgc); %#ok<*NBRAK>
                    myt = (1:length(dgc))*t.dr.derivTime;
                    myt = myt - mean(myt);
                    plot (tx, xc./np,txcrap,xccrap./npcrap,myt,dgc); xlim([-50 50])
                elseif reply==9
                    clf;
                    close all;
                else
                    satisfied2 = true;
                end
            end
                
            
        elseif num == 7
            clf;
            close all;
            disp(' ')
            disp('PLOT COM')
            timeAxis = t0:1:t1;
            [time, centerofmass] = tempset.meanField2vsField1('eti', 'sloc', timeAxis);
            xposition = centerofmass(1,:);
            
            %cleans up invalid data points
            while numel(time) > numel(xposition)
                timeAxis(numel(time))=[];
            end
            
            k=1;
            while k<=numel(time)
                if isnan(xposition(k))
                    centerofmass(:,k)=[];
                    xposition(k)=[];
                    time(k)=[];
                else
                    k=k+1;
                end
            end
            
            dims = get(0,'ScreenSize');
            h=dims(4)-30;
            
            
            fig = figure;
            figs2save(length(figs2save)+1)=fig;
            fignames(length(figs2save))={'composx'};
            delx = xposition - xposition(1,1);
            warning off MATLAB:polyfit:RepeatedPointsOrRescale;
            t=time(1)+25:1:time(numel(time))-25;
            if (time(length(time))-time(1))/3600 < 0.8
                pn = 1;
            else
                pn = 20;
            end
            fit = polyval(polyfit(time,delx,pn),t);
            plot (time,delx,'.',t,fit,'LineWidth',2,'Color','r','MarkerEdgeColor','g');
            xlabel('time');
            ylabel('x pos of CoM');
            axis([time(1) time(length(time)) 0 1]);
            axis 'auto y';
            if pn==1
                legend('data','lin fit','Location','SouthOutside');
            else
                legend('data','poly fit','Location','SouthOutside');
            end
            set(fig,'OuterPosition',[0 h/2+35 500 h/2-5]);
            hold on
            zt = 0:1:t1;
            plot(zt,0,'Color','k');
            hold off;
            
            fig = figure;
            figs2save(length(figs2save)+1)=fig;
            fignames(length(figs2save))={'comvelx'};
            sigma=1;
            xcomvel = deriv(fit,sigma)./deriv(t,sigma);
            plot (t, xcomvel,'LineWidth',2);
            hold on
            plot(t,mean(xcomvel),'LineWidth',2,'Color','r');
            xlabel('time');
            ylabel('x vel of CoM');
            axis([time(1) time(length(time)) 0 1]);
            axis 'auto y';
            if pn==1
                legend('lin fit',['mean vel: ' num2str(mean(xcomvel))],'Location','SouthOutside');
            else
                legend('poly fit',['mean vel: ' num2str(mean(xcomvel))],'Location','SouthOutside');
            end
            set(fig,'OuterPosition',[505 h/2+35 500 h/2-5]);
            plot(zt,0,'Color','k');
            hold off;
            
            fig = figure;
            figs2save(length(figs2save)+1)=fig;
            fignames(length(figs2save))={'composy'};
            yposition = centerofmass(2,:);
            dely = yposition - yposition(1,1);
            warning off MATLAB:polyfit:RepeatedPointsOrRescale;
            if (time(length(time))-time(1))/3600 < 0.8
                pn = 1;
            else
                pn = 20;
            end
            fit = polyval(polyfit(time,dely,pn),t);
            plot (time, dely,'.',t,fit,'LineWidth',2,'Color','r','MarkerEdgeColor','g');
            xlabel('time');
            ylabel('y pos of CoM');
            axis([time(1) time(length(time)) 0 1]);
            axis 'auto y';
            if pn==1
                legend('data','lin fit','Location','SouthOutside');
            else
                legend('data','poly fit','Location','SouthOutside');
            end
            set(fig,'OuterPosition',[0 30 500 h/2-5]);
            hold on
            plot(zt,0,'Color','k');
            hold off;
            
            fig = figure;
            figs2save(length(figs2save)+1)=fig;
            fignames(length(figs2save))={'comvely'};
            ycomvel = deriv(fit, sigma)./deriv(t, sigma);
            plot (t, ycomvel,t,mean(ycomvel),'LineWidth',2);
            hold on
            plot(t,mean(ycomvel),'LineWidth',2,'Color','r');
            xlabel('time');
            ylabel('y vel of CoM');
            axis([time(1) time(length(time)) 0 1]);
            axis 'auto y';
            if pn==1
                legend('lin fit',['mean vel: ' num2str(mean(ycomvel))],'Location','SouthOutside');
            else
                legend('poly fit',['mean vel: ' num2str(mean(ycomvel))],'Location','SouthOutside');
            end
            set(fig,'OuterPosition',[505 30 500 h/2-5]);
            plot(zt,0,'Color','k');
            hold off;

        elseif num == 8
            filenames='';
            for k=1:numel(tempset.expt)
                filenames = [filenames ' ' tempset.expt(k).fname(ind+1:numel(tempset.expt(k).fname))]; %#ok<*AGROW>
            end
            disp(' ')
            disp('  Experiment Info')
            disp(['  |     file(s):' filenames]);
            disp(['  |    # tracks: ' num2str(numel([tempset.expt.track]))]);
            if ~multi
                disp(['  | time elapse: ' num2str(t1-t0) ' sec (' num2str(t0) '~' num2str(t1) ' sec)']);
                disp(['  |   dimension: ' num2str(x1-x0) ' x ' num2str(y1-y0) ' px']);
                disp(['  |   queue pos: ' num2str(i) ' of ' num2str(numel(exs.expt)) ' experiments']);
            else
                disp('  |   GROUP ANALYSIS IN PROGRESS');
            end
            
        elseif num==9
            clf;
            close all;
            disp(' ')
            disp('CLEAR/CLOSE ''EM')
            disp(' ')
            
        elseif num==-1
            if ~multi
                disp(' ')
                rex = ExperimentSet.fromFiles(ex.fname,'minpts',200);
                exs.expt(i) = rex.expt(1);
                ex = exs.expt(i);
                tempset.expt=[ex];
                segmented = false;
                disp(' ')
                disp('EXPERIMENT RELOADED.')
            else
                disp(' ')
                disp('Reloading not allowed for group analysis!')
                disp(' ')
            end
        end
        
        if ~isempty(figs2save)
            reply = input('\n  Save figures?\n  | 0. No\n  | 1. Yes\n  | Your Response: ');
            if reply == 0
                clf;
                close all;
            else
                disp(' ')
                for q=1:length(figs2save)
                    figfile=[figsdir '\' char(fignames(q)) '1'];
                    copyn = 1;
                    while exist([figfile '.fig'])>0 %#ok<*EXIST>
                        copyn = copyn+1;
                        figfile = [figfile(1:length(figfile)-1) num2str(copyn)];
                    end
                    saveas(figs2save(q),[figfile '.fig']);
                    disp(['  | ' char(fignames(q)) num2str(copyn) '.fig saved!'])
                end
            end
        end
    end
end

disp(' ')
disp(line)
disp(' ')

disp('We''re done! WOOT. :D')