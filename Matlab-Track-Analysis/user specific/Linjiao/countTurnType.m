% count turn type

a=cell(1,100);
b=[];
num=0; % num of turn types
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        for k=1:length(eset.expt(i).track(j).reorientation)
            reo=eset.expt(i).track(j).reorientation(k);
            reoseq=reo.turnsequence;
            new=1; %a new turn type
            if num>0
               for n=1:num
                   if isequal(a{1,n},reoseq);
                       b(n)=b(n)+1;
                       new=0;
                   end
               end
            end
            if new
                num=num+1;
                a{1,num}=reoseq;
                b(num)=1;
            end
        end
    end
end

[bb, ind]=sort(b,'descend');
% for i=1:length(bb)
%     cc(i)=num2str(a{1,ind(i)});
% end
bar(bb(1:10));

           