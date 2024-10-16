fignum=figure(); 
clf(fignum);  %create a figure and clear it
hold on;
for i=1:numel(eset)
eset(i).executeTrackFunction('plotPath', 'displacement'); % plot all paths
plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
axis([-1200 1200 -1000 1000]);
set(gca, 'XTick',[-900,-450,0,450,900]);
set(gca,'ycolor',[1 1 1]) 
set(gca,'ytick',[]);
end
hold off;