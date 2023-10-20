function gq = vocppmFromSetPoint (gqSetPoint, amplitude, sigmaInSec, offsetInSec)
%function gq = vocppmFromSetPoint (gqSetPoint, amplitude, sigmaInSec, offsetInSec)
%I think amplitude should be ppm on meter / setpoint
%
existsAndDefault('sigmaInSec', 33.7563);
existsAndDefault('offsetInSec', 8.6659);

et = gqSetPoint.xData.et;
dt = median(diff(et));
sigma = sigmaInSec/dt;
offset = offsetInSec/dt;

kernel = gaussWithOffset(amplitude, sigma, offset);
gq = gqSetPoint;
gq.yData = conv(gq.yData, kernel, 'same');
gq.fieldname = 'vocppm';
