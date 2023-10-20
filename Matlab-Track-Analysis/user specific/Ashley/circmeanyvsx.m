function [x,meany,standarderror,standarddeviation] = circmeanyvsx (xdata, ydata, xaxis)
%function [x,meany,standarderror,standarddeviation] = circmeanyvsx (xdata, ydata, xaxis)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each interval in xaxis (xdata >= xaxis(j), <= xaxis(j+1)), 
%meany is the mean of all ydata with xdata in
%that interval
%
%standard error is the standard deviation of the data in the bin divided by
%the square root of the number of elements

for j = 1:(length(xaxis) - 1)
    inds = find(xdata >= xaxis(j) & xdata < xaxis(j+1));
    x(j) = circ_mean(xdata(inds));
    meany(j) = circ_mean(ydata(inds));
    standarderror(j) = circ_std(ydata(inds))/sqrt(length(inds));
    standarddeviation(j) = circ_std(ydata(inds));
end
