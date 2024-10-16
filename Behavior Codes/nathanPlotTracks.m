trackCount = 0;
startX = 0;
avgStart = 0;
startTimeLimit = 1000;

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
        for k = 1:length(eset(i).expt(j).track);
        if eset(i).expt(j).track(1,k).startFrame < startTimeLimit
            trackCount = trackCount + 1;
            startX = startX + eset(i).expt(j).track(1,k).pt(1,1).loc(1);
        end
        end
    end
end

avgStart = startX/trackCount;

fignum=figure();
clf(fignum);  %create a figure and clear it
hold on;
eset.executeTrackFunction('plotPath'); % plot all paths
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
axis([0 2500 0 2000]);
yL = get(gca,'YLim');
line([avgStart avgStart],yL,'Color','r');
hold off;