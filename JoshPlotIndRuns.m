

for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
        % for k = 1:length(eset(i).expt(j).track);
            k =1; 
            fignum=figure();
            clf(fignum);  %create a figure and clear it
            
            for l = 1:length(eset(i).expt(j).track(k).runs);
                
                eset(i).expt(j).track(1,k).run(l).plotPath();
                % plot (0,0, 'r.', 'MarkerSize', 20); % put a red dot at 0,0
                axis([0 2500 0 2000]);
                yL = get(gca,'YLim');
                avgStart=eset(i).expt(j).track(1,k).pt(1,1).loc(1);
                line([avgStart avgStart],yL,'Color','r');
            end
        % end
    end
end

