len=length(eset.expt);
while length((eset.expt(len).track))==0
      len=len-1;
end
len1=len;
for i=2:len1
    if length(eset.expt(len1-i+1).track)==0
        for j=(len1-i+1):len-1
            eset.expt(j)=eset.expt(j+1);
        end
        len=len-1;
    end
end
eset.expt=eset.expt(1:len);

        