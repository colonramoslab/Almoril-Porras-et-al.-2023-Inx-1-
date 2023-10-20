function [fignum]=ExpPlotTracks2(expt);
% Faster than old nathanPlotTracks,
% Can pick color by using colChoice or use the angle of migration to color (slower)

trackCount = 0;
startX = 0;
avgStart = 0;
startTimeLimit = 1000;
angleColor=0; % 1 to color by theta angle of run,
spotSize=1;
colChoice=rand([1,3]);

%varargin = assignApplicable, blah, blah blah


fignum=figure();
clf(fignum);  %create a figure and clear it
hold on;
plot (0,0, 'r.', 'MarkerSize', 5); % put a red dot at 0,0
axis([0 2500 0 2000]);
yL = get(gca,'YLim');
if length(expt.track)>1
    if angleColor==1
        for j = 1:length(expt);
            for k = 1:length(expt(j).track);
                trk=expt(j).track(k);
                x=trk.dq.sloc(1,:); y=trk.dq.sloc(2,:); z=trk.dq.theta(:);
                scatter(x, y, spotSize, -abs(z),'filled');
                colormap(jet);
            end
        end
    elseif ~isempty(expt.track(1).dq)
        for j = 1:length(expt);
            for k = 80;
                trk=expt(j).track(k);
                x=trk.dq.sloc(1,:); y=trk.dq.sloc(2,:);
                plot(x, y); %,'color',colChoice);
                colormap(jet);
            end
        end
    else
        for j = 1:length(expt);
            for k=1:length(expt.track)
                plotPath(expt(j).track(k));
            end
        end
    end
    
    
    
    avgStart = startX/trackCount;
    line([avgStart avgStart],yL,'Color','r');
    hold off;
end
end