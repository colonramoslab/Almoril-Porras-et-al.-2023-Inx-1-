function eset_statistics = calculateStatisticsOfEset(eset, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (length(eset) > 1)
    for j = 1:length(eset)
        es(j) = calculateStatisticsOfEset(eset(j), varargin{:});
    end
    eset_statistics = es;
    return;
end

timerange = [];
directions = [0 180 -90 90];
varargin = assignApplicable(varargin);

es.directions = directions;

es.runStartTime = eset.gatherFromSubField ('run', 'eti', 'position', 'start');
es.runEndTime = eset.gatherFromSubField ('run', 'eti', 'position', 'end');

es.reoStartTime = eset.gatherFromSubField ('reorientation', 'eti', 'position', 'start');
%es.reoEndTimes = eset.gatherFromSubField ('reorientation', 'eti', 'position', 'end');

es.runStartTheta = eset.gatherSubField ('run', 'startTheta');
es.runEndTheta = eset.gatherSubField ('run', 'endTheta');
es.runMeanTheta = eset.gatherSubField('run' ,'meanTheta');

es.reoPrevDir = eset.gatherSubField ('reorientation', 'prevDir');
es.reoNextDir = eset.gatherSubField ('reorientation', 'nextDir');
es.reoNumHS = eset.gatherSubField ('reorientation', 'numHS');

es.headSwingTime = eset.gatherFromSubField ('headSwing', 'eti', 'position', 'start');
es.firstHeadSwingTime = eset.gatherFromSubField ('firsths', 'eti', 'position', 'start');
es.lastHeadSwingTime = eset.gatherFromSubField ('lasths', 'eti', 'position', 'start');

es.headSwingAccepted = eset.gatherSubField ('headSwing', 'accepted');
es.headSwingPrevDir = eset.gatherSubField ('headSwing', 'prevDir');
es.firstHeadSwingPrevDir = eset.gatherSubField ('firsths', 'prevDir');
es.firstHeadSwingAccepted = eset.gatherSubField ('firsths', 'accepted');
es.lastHeadSwingAccepted = eset.gatherSubField ('lasths', 'accepted');
es.headSwingValid = eset.gatherSubField ('headSwing', 'valid');
es.firstHeadSwingValid = eset.gatherSubField ('firsths', 'valid');
es.lastHeadSwingValid = eset.gatherSubField ('lasths', 'valid');

eti = eset.gatherField('eti');

if (isempty(timerange))
    timerange = [min(eti)-1 max(eti)+1];
end
es.timerange = timerange;

it = eset.gatherSubField('dr','interpTime');
dt = median(it);
if (any(it ~= dt))
    disp (['warning:  eset does not have homogenous interpolation times, instead range from ' num2str(min(it)) ' to ' num2str(max(it))]);
end

es.maxNumAnimals = zeros(size(eset.expt));
es.maxNumAnimalsInTimeWindow = es.maxNumAnimals;
for j = 1:length(eset.expt)
    et = eset.expt(j).gatherField('eti');
    es.maxNumAnimals(j) = max(histc(et, min(et):1:max(et)))*dt;
    es.maxNumAnimalsInTimeWindow(j) = max(histc(et, min(timerange):1:max(timerange)))*dt;
end
es.numExpts = length(eset.expt);
es.numAnimals = sum(es.maxNumAnimalsInTimeWindow);
es.animalTime = sum(eti >= min(timerange) & eti <= max(timerange))*dt;

for j = 1:length(directions)
    es.numRunsInDirection(j) = nnz(cos(es.runMeanTheta - deg2rad(directions(j))) > 1/sqrt(2) & es.runStartTime >= min(timerange) & es.runStartTime <= max(timerange));
    es.numReosFromDirection(j) = nnz(cos(es.reoPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.reoStartTime >= min(timerange) & es.reoStartTime <= max(timerange));
    es.numReosWithHSFromDirection(j) = nnz(cos(es.reoPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.reoStartTime >= min(timerange) & es.reoStartTime <= max(timerange) & es.reoNumHS > 0);
    es.numHSFromDirection(j) = nnz(cos(es.headSwingPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.headSwingTime >= min(timerange) & es.headSwingTime <= max(timerange) & es.headSwingValid);
end

eset_statistics = es;

