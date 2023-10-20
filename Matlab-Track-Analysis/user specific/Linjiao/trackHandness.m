handness=[]; % find the handness of the tracks based all the single omega turns
numberO=[];
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        left=0;
        right=0;
        for k=1:length(eset.expt(i).track(j).reorientation)
            reo=eset.expt(i).track(j).reorientation(k);
            if isequal(reo.turnsequence,[-1]);
                deltatheta=reo.thetaOut-reo.thetaIn;
                deltatheta=rad2deg(deltatheta);
                if deltatheta>180
                    deltatheta=deltatheta-360;
                elseif deltatheta<-180
                    deltatheta=deltatheta+360;
                end
                if deltatheta>0
                    left=left+1;
                elseif deltatheta<0
                    right=right+1;
                end
            end
        end
        hand=left/(left+right);
        handness=[handness, hand];
        numberO=[numberO, left+right];
    end
end
                
                
            