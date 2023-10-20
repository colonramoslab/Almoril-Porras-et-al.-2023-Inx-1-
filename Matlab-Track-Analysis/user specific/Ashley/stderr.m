function [std_err]= stderr(array)

std_err = std(array)/sqrt(length(array));

end