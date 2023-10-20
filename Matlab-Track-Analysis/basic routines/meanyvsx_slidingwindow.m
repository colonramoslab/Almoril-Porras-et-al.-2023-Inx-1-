function [x,meany,standarderror,standarddeviation,sumy, numy] = meanyvsx_slidingwindow (xdata, ydata, xcenter, bindim, windowType, polar)
%function [x,meany,standarderror,standarddeviation, sumy, numy] = meanyvsx (xdata, ydata, xaxis, xcenter, bindim, windowType, polar)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each point in xcenter, we take the weighted average of ydata with
%weighting determined by windowType and bindim
%
%if windowType is 'gaussian' (default), the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) [bindim = FWHM]

%if windowType is 'step', the weighting w_j is 1 if abs(x_j - xc) < bindim/2
%and 0 otherwise [bindim = bin width]
%
%if windowType is 'halfgaussian', the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) * (sign(x_j - xc) == sign(bindim)) :  in other words, if
%bindim > 0 we only consider xj > xc, and if bindim is < 0 we only consider
%xj < xc
%
%if polar is true, we add or sutbract multiples of 2*pi to the xdata to lie
%within +/-pi of the given center
%meany is the mean of all ydata with xdata in
%that interval
%
%
%standard error is the standard deviation of the data in the bin divided by
%the square root of the number of elements
existsAndDefault('polar', false);
sumy = zeros([size(ydata,1) length(xcenter)]);
meany = sumy;
numy = zeros(1, length(xcenter));
x = numy;
standarddeviation = meany;
standarderror = meany;

switch(lower(windowType))
    case 'gaussian'
        weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2);
    case 'halfgaussian'
        weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2).*(sign(xj-xc) == sign(bindim));
    case 'step'
        weightFun = @(xd,xc) heaviside(bindim/2 - abs(xd-xc));
    otherwise
        disp('windowType must be gaussian or step');
        return;
        
end
valid = isfinite(xdata) & all(isfinite(ydata),1);
xdata = xdata(valid);
ydata = ydata(:,valid);

for j = 1:length(xcenter)
    if (polar)
        xdata = mod(xdata - xcenter(j) + pi, 2*pi) - pi + xcenter(j);
    end
    
    w = weightFun(xdata,xcenter(j));
   
    numy(j) = sum(w);
    x(j) = sum(w.*xdata)/numy(j);
    w = repmat(w, size(ydata, 1), 1);
    
    sumy(:,j) = sum(w.*ydata, 2);
    meany(:,j) = sumy(j) / numy(j);
    standarddeviation(:,j) = sqrt(sum(w.*(ydata - repmat(meany(:,j), 1, size(ydata, 2))).^2)/numy(j));
    standarderror(:,j) =  standarddeviation(:,j)/sqrt(numy(j));
end


