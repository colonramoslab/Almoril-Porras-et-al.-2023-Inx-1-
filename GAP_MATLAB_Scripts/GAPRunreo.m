function []= GAPRunreo(dDir)

fileL=dir([dDir, '/*.mat']);


for i=1:length(fileL)
% Load & segment file
load(fullfile(fileL(i).folder,fileL(i).name));
segmentTracks(experiment_1);
expt=experiment_1;
segmentOptions=WormSegmentOptions;
expt(1).segmentTracks(segmentOptions);
for j = 1:length(expt)
            minXPosition = expt(j).track(1,1).dq.sloc(1,1);
            minYPosition = expt(j).track(1,1).dq.sloc(2,1);
            maxXPosition = expt(j).track(1,1).dq.sloc(1,1);
            maxYPosition = expt(j).track(1,1).dq.sloc(2,1);
            for k = 1:length(expt(j).track)
                for l = 1:length(expt(j).track(1,k).dq.sloc)
                    if expt(j).track(1,k).dq.sloc(1,l)<minXPosition
                        minXPosition = expt(j).track(1,k).dq.sloc(1,l);
                    end
                    if expt(j).track(1,k).dq.sloc(2,l)<minYPosition
                        minYPosition = expt(j).track(1,k).dq.sloc(2,l);
                    end
                    if expt(j).track(1,k).dq.sloc(1,l)>maxXPosition
                        maxXPosition = expt(j).track(1,k).dq.sloc(1,l);
                    end
                    if expt(j).track(1,k).dq.sloc(2,l)>maxYPosition
                        maxYPosition = expt(j).track(1,k).dq.sloc(2,l);
                    end
                end
            end
end
fignum=figure();
clf(fignum);  %create a figure and clear it
hold on;
axis([minXPosition maxXPosition minYPosition maxYPosition]);
yL = get(gca,'YLim');
for j = 1:length(expt)
            for k = 1:length(expt(j).track)
                trk=expt(j).track(k);
                x=trk.dq.sloc(1,:); y=trk.dq.sloc(2,:);
                plot(x, y,'k');
            end
end

for j = 1:length(expt)
            for k = 1:length(expt(j).track)
                for l = 1:length(expt(j).track(1,k).run)
                    colorchange = rem(l,2);
                    runStart=expt(j).track(1,k).run(1,l).startInd;
                    runEnd=expt(j).track(1,k).run(1,l).endInd;
                    trk=expt(j).track(k);
                    x=trk.dq.sloc(1,runStart:runEnd); y=trk.dq.sloc(2,runStart:runEnd);
                        plot(x, y,'color',[(0.35+(0.35*colorchange)) (0.35+(0.35*colorchange)) (0.35+(0.35*colorchange))]);
                end
            end
end
for j = 1:length(expt)
            for k = 1:length(expt(j).track)
                for l = 1:length(expt(j).track(1,k).reorientation)
                    reoStart=expt(j).track(1,k).reorientation(1,l).startInd;
                    reoEnd=expt(j).track(1,k).reorientation(1,l).endInd;
                    trk=expt(j).track(k);
                    x=trk.dq.sloc(1,reoStart:reoEnd); y=trk.dq.sloc(2,reoStart:reoEnd);
                    plot(x, y,'r');
                end
            end
end
for j = 1:length(expt)
            for k = 1:length(expt(j).track)
                for l = 1:length(expt(j).track(1,k).reorientation)
                    for m = 1:length(expt(j).track(1,k).reorientation(1,l).sharpTurn)
                        turnStart=expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).startInd;
                        turnEnd=expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).endInd;
                        trk=expt(j).track(k);
                        x=trk.dq.sloc(1,turnStart:turnEnd); y=trk.dq.sloc(2,turnStart:turnEnd);
                        if expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).typeCode == -1 % -1 = 'omega turn'; 0 = 'double reverse or blip'; 1 = 'reversal'; 2 = 'second reversal'
                                plot(x, y,'g');
                        elseif expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).typeCode == 0
                                plot(x, y,'y');
                        elseif expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).typeCode == 1
                                plot(x, y,'b');
                        elseif expt(j).track(1,k).reorientation(1,l).sharpTurn(1,m).typeCode == 2
                                plot(x, y,'c');
                        end
                    end
                end
            end
end
            set(fignum,'Units','Inches');
            pos=get(fignum,'Position');
            set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
                'PaperSize',[pos(3),pos(4)]);
            outName=strcat(fullfile(fileL(i).folder,fileL(i).name),'_run_vs_reo_reversalCovThresh_1_6.pdf');
            print(fignum,outName,'-dpdf','-r150','-painters');
            close(fignum);
            minXPosition = 0;
            minYPosition = 0;
            maxXPosition = 0;
            maxYPosition = 0;
            trk = 0;
end
end