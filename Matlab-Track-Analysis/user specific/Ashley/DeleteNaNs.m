function [newarray] = DeleteNaNs (arraywNaN)
%function [newarray] = DeleteNaNs (arraywNaN)
%
%this function finds the mean value of the array 'arraywNaN' and outputs it

newarray=[];
count=0;
for n=1:length(arraywNaN)
    if(arraywNaN(n)<Inf)
        count=count+1;
        newarray(count)=arraywNaN(n);
    end
end



end