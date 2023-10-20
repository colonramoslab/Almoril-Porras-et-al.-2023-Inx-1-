function [lPass, contLengths] = stretchCounter(passVals)
%UNTITLED Create vector with each logical value replaced by length of
%continuity of this logical value in array. countLengths holds the lengths
%for each continuous run of logical trues.

% first, encode length of a given stretch

lPass=NaN(size(passVals));
contLengths=[];
holdVal=0;
for ii=1:length(passVals)
    if passVals(ii)==0
        % assign hole value
        lPass(ii)=0;
        % assign previous stretch based on accumulated 'holdVal'
        if holdVal>0
            lPass(ii-holdVal:ii-1)=holdVal;
              %reset holdVal
            contLengths=[contLengths,holdVal];
            holdVal=0;
        end
        
    else
        holdVal=holdVal+1;
    end
end
end

