function [cm, ccov, cmr] = meanAndCovEMStep(sc, point_source, postprob, numclusters)
%finds the mean and covariance for clustered data, based on posterior
%probability
%function [cm, ccov] = meanAndCovEMStep(sc, point_source, postprob)
%
%notation: K is number of data fields in sc
%outputs:
%   CM: Kx1 vector that represents the center of the cluster
%           if sc is an array, then CM is a cell of vectors
%   CCOV: KxK covariance matrix; if sc is an array then CCOV is a cell
%inputs:
%   SC < SegmentationCluster
%   POINT_SOURCE < Track, Experiment, ExperimentSet
%   POSTPROB: posterior probability given by posteriorProbabilities(scSet,
%           track);
%           if sc is a single cluster, postprob must be passed
%           if sc is an array and postprob is not passed, then we calculate
%           posterior probability using point_source and sc

if (length(sc) > 1)
    if (~exist('postprob', 'var') || isempty(postprob))
        postprob = sc.posteriorProbabilities(point_source);
    end
    for j = 1:length(sc)
        [out1,out2] = meanAndCovEMStep(sc(j), point_source, postprob(j,:));
        cm{j} = out1; %#ok<AGROW>
        ccov{j} = out2; %#ok<AGROW>
    end
    return;
end

existsAndDefault('numclusters', 1);


if (isa(point_source, 'Track'))
    data = zeros(length(sc.datafields), length(point_source.getDerivedQuantity('eti')));

    for j = 1:length(sc.datafields)
        op = sc.operation{j};
        data(j,:) = op(point_source.getDerivedQuantity(sc.datafields{j}));
    end
else
    if (ismethod(point_source, 'gatherField'))
        data = zeros(length(sc.datafields), length(point_source.gatherField('eti')));

        for j = 1:length(sc.datafields)
            op = sc.operation{j};
            data(j,:) = op(point_source.gatherField(sc.datafields{j}));
        end
    else
        warning('sc:mac', 'point source must be a track, experiment, or experiment set');
        cm = [];
        ccov = [];
        return;
    end
end

pp = repmat(postprob, [size(data,1), 1]);
gmfit = gmdistribution.fit((pp.*data)', numclusters,'Regularize',1E-6);
cmr = gmfit.PComponents;
ccov = gmfit.Sigma;
cm = gmfit.mu;

%{
cm = mean(pp.*data, 2);
data = data - repmat(cm, [1, size(data,2)]);
ccov = cov(transpose(pp.*data));
%}
