function [new_ang] = ClosestAng(angin)

test=[-pi -pi/2 0 pi/2 pi];
if(length(angin)==1)
    [minang, minint]=min(abs(test-angin));
    new_ang=test(minint);
    if(new_ang==pi)
        new_ang=-pi;
    end
else
    new_ang(1:length(angin))=0;
    for j=1:length(angin)
        [minang, minint]=min(abs(test-angin(j)));
        new_ang(j)=test(minint);
    end
    new_ang(new_ang==pi)=-pi;
end

end