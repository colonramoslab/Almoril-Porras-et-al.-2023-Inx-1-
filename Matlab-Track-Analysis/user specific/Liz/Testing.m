count=0;
for j=1:length(onTimeCount)-1
    if onTimeCount(j)==1
        onTimeCount(j+1)=onTimeCount(j)+1;
    end
end
