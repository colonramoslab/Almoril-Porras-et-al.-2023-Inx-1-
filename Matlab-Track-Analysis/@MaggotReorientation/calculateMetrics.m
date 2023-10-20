function calculateMetrics(reo, prevRun, nextRun)
%function calculateMetrics(reo)

if (~isempty(prevRun) && isa(prevRun, 'Run'))
    startInd = prevRun.endInd + 1;
else
    startInd = [];
end
if (~isempty(nextRun) && isa(nextRun, 'Run'))
    endInd = nextRun.startInd - 1;
else
    startInd = [];
end


reo.numHS = length(reo.headSwing);
if (~isempty(reo.headSwing))
    reo.startInd = min([reo.headSwing.startInd]);
    reo.endInd = max([reo.headSwing.endInd]);
else
    reo.startInd = startInd;
    reo.endInd = endInd;
end
reo.inds = reo.startInd:reo.endInd;

if (~isempty(reo.prevRun))
    reo.prevDir = reo.prevRun.endTheta;
end
if (~isempty(reo.nextRun))
    reo.nextDir = reo.nextRun.startTheta;
end

for j = 1:length(reo.headSwing)
    reo.headSwing(j).prevDir = reo.prevDir;
    reo.headSwing(j).nextDir = reo.nextDir;
end
