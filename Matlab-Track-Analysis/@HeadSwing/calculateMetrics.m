function calculateMetrics(hs)

    hs.track.calculateDerivedQuantity({'sbodytheta', 'shead', 'stail', 'smid'});
    last = find(hs.track.isrun(hs.startInd:hs.endInd),1, 'first');
    if (~isempty(last))
        hs.endInd = hs.startInd - 1 + last;
    end
        
    hs.inds = hs.startInd:hs.endInd;
    btheta = hs.track.dq.sbodytheta(hs.inds);
    mh = hs.track.dq.shead(:,hs.inds) - hs.track.dq.smid(:,hs.inds);
    tm = hs.track.dq.smid(:,hs.inds) - hs.track.dq.stail(:,hs.inds);
    
    
    [blah,I] = max(abs(btheta));
    hs.maxInd = hs.inds(I);
    hs.maxTheta = btheta(I);
    hs.sign = sign(hs.maxTheta);
    hs.headDir = atan2(mh(2,I), mh(1,I));
    hs.tailDir = atan2(tm(2,I), tm(1,I));
    hs.accepted = hs.track.isrun(hs.endInd);
    hs.valid = all(hs.track.getDerivedQuantity('ihtValid', false, 'inds', hs.inds));
    if (~isempty(hs.nextRun))
        hs.nextDir = hs.nextRun.startTheta;
    end
    if (~isempty(hs.prevRun))
        hs.prevDir = hs.prevRun.endTheta;
    end
    %{
    hs.prevRun = hs.track.run(find([hs.track.run.stop] < hs.endInd, 1, 'last'));
    hs.nextRun = hs.track.run(find([hs.track.run.start] < hs.startInd, 1, 'first'));
    if (~isempty(hs.prevRun))
        hs.prevDir = hs.prevRun.endTheta;
    end    
    hs.nextDir
    %}
end