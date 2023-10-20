function [eset] = DeleteTracksThatAreSlow(eset)
%function [eset] = DeleteTracksThatAreSlow(eset)
%
%looks at the distribution of avg speeds for tracks in each experiment and
%deletes those that are too slow

expt=ExperimentSet();

for j=1:length([eset.expt])
    expt=eset.expt(j);
    [s]=expt.gatherField('speed','mean');
    outlier = s-mean(s) < -2*std(s);
    expt.track=expt.track(find(~outlier));
    eset.expt(j)=expt;
    
end


end

