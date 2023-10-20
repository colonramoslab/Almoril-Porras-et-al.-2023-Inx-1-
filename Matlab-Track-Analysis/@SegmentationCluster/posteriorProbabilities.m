function postprob = posteriorProbabilities(scSet, point_source)
%the posterior probability that each ind in POINT_SOURCE is part of each cluster
%function postprob = posteriorProbabilities(scSet, POINT_SOURCE)
%
%POINT_SOURCE is one of (Track, Experiment, ExperimentSet) 
%likelihood(j,ind) =  scSet(j).priorProb * scSet(j).likelihood(point_source)
%postrob(j,ind) = likelihood(j,ind)/sum(likelihood(:,ind));

if (isa(point_source, 'Track'))
    dim2 = length(point_source.getDerivedQuantity('eti'));
else
    if (ismethod(point_source, 'gatherField'))
        dim2 = length(point_source.gatherField('eti'));
    else
        postprob = [];
        return;
    end
end
lk = zeros(length(scSet), dim2);
for j = 1:length(scSet)
    lk(j,:) = scSet(j).likelihood(point_source)*scSet(j).priorProbability;
end
nf = repmat(sum(lk,1),[size(lk,1) 1]);
postprob = lk./nf;
