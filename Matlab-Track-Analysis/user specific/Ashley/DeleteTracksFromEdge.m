function [eset] = DeleteTracksFromEdge(eset, limits)
%function [eset] = DeleteTracksFromEdge(eset, limits)
%
%Delete all of the tracks that start from the edge, edge is defined by lims
%limits = [xstart xend ystart yend]

% x_start_lim=400;
% x_end_lim=1950;
% y_start_lim=150;
% y_end_lim=1750;

for n=1:length([eset.expt])
    expt=eset.expt(n);
    count=0;
    for j=1:length([expt.track])
        count=count+1;
        if((expt.track(count).pt(1).loc(1)<limits(1))||(expt.track(count).pt(1).loc(1)>limits(2)))
            expt.track(count)=[];
            count=count-1;
        elseif((expt.track(count).pt(1).loc(2)<limits(3))||(expt.track(count).pt(1).loc(2)>limits(4)))
            expt.track(count)=[];
            count=count-1;
        end
    end
    
    eset.expt(n)=expt;
end


end

