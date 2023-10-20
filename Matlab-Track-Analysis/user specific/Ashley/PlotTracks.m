function [] = PlotTracks(eset)
%function [] = PlotTracks(eset)
%
%plots all of the tracks for each experiment on 1 graph with the background

back=imread('C:\Users\Ashley\Documents\Data\2010\Mar\19\0930\0930_background\0930_background1.jpg');

for n=1:length(eset)
    expt=eset(n).expt;
    figure(n);
    imagesc(back);
    colormap('gray');
    hold on;
    expt.executeTrackFunction('plotPath');
    hold off;
    
end

end