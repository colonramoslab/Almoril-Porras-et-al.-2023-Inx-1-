function lk = likelihood (wsc, point_source)
%P(data|cluster)
%function lk = likelihood (sc, point_source)
%
%outputs: 
%   lk - for each interpolated point in point_source, the likelihood of getting
%        that data point from the cluster
%inputs:
%   sc < SegmentationCluster
%   point_source < Track, Experiment, or ExperimentSet

if (isempty(wsc.knownPoints) || isempty(wsc.knownPoints(1).track) || isempty(wsc.knownPoints(1).inds))
    lk = 0;
    return;
end

if isempty(wsc.operation)
    wsc.operation = {@(x) x};
    for j = 1:length(wsc.datafields)
        wsc.operation{j} = @(x) x;
    end
end

if (isa(point_source, 'Track'))
    data = zeros(length(wsc.datafields), length(point_source.getDerivedQuantity('eti')));

    for j = 1:length(wsc.datafields)
        op = wsc.operation{j};
        data(j,:) = op(point_source.getDerivedQuantity(wsc.datafields{j})) - wsc.clustMean(j);
    end
else
    if (ismethod(point_source, 'gatherField'))
        data = zeros(length(wsc.datafields), length(point_source.gatherField('eti')));

        for j = 1:length(wsc.datafields)
            op = wsc.operation{j};
            data(j,:) = op(point_source.gatherField(wsc.datafields{j})) - wsc.clustMean(j);
        end
    else
        lk = [];
        return;
    end
end
prefactor = 1/ ((2 * pi)^(size(data,1)/2) * det(wsc.clustCov));
ic = inv(wsc.clustCov);

lk = prefactor * exp(-0.5*dot(data, ic * data));
