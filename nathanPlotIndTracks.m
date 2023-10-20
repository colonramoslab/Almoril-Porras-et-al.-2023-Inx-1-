

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
        for k = 1:length(eset(i).expt(j).track);
            fignum=figure();
            clf(fignum);  %create a figure and clear it
            eset(i).expt(j).track(1,k).plotPath();
           % plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
            axis([0 2500 0 2000]);
            yL = get(gca,'YLim');
            avgStart=eset(i).expt(j).track(1,k).pt(1,1).loc(1);
            line([avgStart avgStart],yL,'Color','r');
        end
    end
end

