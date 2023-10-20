function setParams(s, thetaStep, biasDirection, baseRate, rateDelta, baseMag, magDelta, reoBias)
%BiasedRandomWalkSimulator/
%function setParams(s, thetaStep, biasDirection, baseRate, rateDelta,
%baseMag, magDelta, reoBias);
%
%set parameters for Biased Random Walk
%
%thetaStep: how fine grained the step size should be (in radians)
%biasDirection: the direction to which the random walk will be biased
%baseRate: the base rate (in inverse seconds) of reorientations
%rateDelta: the fractional change when heading in the wrong/right direction
%           rate = baseRate * (1 - rateDelta*cos(dir - biasDirection));
%baseMag: the base standard deviation for the size of a turn (in radians)
%magDelta: the fractional change when heading in the wrong/right direction
%          sd = baseMag * (1 - magDelta*cos(dir - biasDirection));
%ctrBias: the amount to which the reorientation is biased 
%       0 = no bias; 1 = always turns in the right direction

s.biasDirection = biasDirection;
s.thetaAxis = (-thetaStep):thetaStep:(2*pi + thetaStep);
s.reoRate = baseRate * (1 - rateDelta*cos(s.thetaAxis - biasDirection));
s.reoMag = baseMag * (1 - magDelta*cos(s.thetaAxis - biasDirection));
s.reoBias = repmat(reoBias, size(s.thetaAxis)); 
