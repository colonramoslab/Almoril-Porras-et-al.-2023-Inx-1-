function sm = fromSegmentedTracks(src, varargin)
%converts a set of segmented tracks into a fully labeled segmentation model
%
%function sm = fromSegmentedTracks(fn, src, varargin)
%src is either ExperimentSet(s), Experiment(s), or Track(s)

    if (isa(src, 'ExperimentSet'))
        sm = SegmentationModel.fromSegmentedTracks([src.expt], varargin{:});
        sm.eset = src;
        return;
    end
    if (isa(src, 'Experiment'))
        sm = SegmentationModel.fromSegmentedTracks([src.track], varargin{:});
        return;
    end
    if (~isa(src, 'Track'))
        error ('I need either an experiment(s), experiment set, or track(s) to work with here');
        return;
    end

    sm = SegmentationModel.twelveStateWormModel;
    for j = 1:length(src)
        addSegmentedTrackToModel(sm, src(j));
    end
end

function addSegmentedTrackToModel (sm, track)

    st = track.sharpTurn;
    inreverse = false;
    maxBackwardRun = 15 / track.dr.interpTime;
    for j = 1:length(st)
    
        switch (st(j).typeCode)
            case -1 %omega turn
                if (inreverse && st(j).startInd - st(j-1).endInd < maxBackwardRun)
                    addOmegaTurn(sm.segmentationClusters, st(j), st(j-1));
                else
                    addOmegaTurn(sm.segmentationClusters, st(j));
                end
                inreverse = false;        
            case 0 
                inreverse = false;
                %blip
            %{
            if checkForOmegaTurn (st(j))
                addOmegaTurn(sm.segmentationClusters, st(j));
                inreverse = false;
            end
            %}
            
            case 1 %reversal
                
                %check time since last to avoid seriously long backwards runs
                if (inreverse && st(j).startInd - st(j-1).endInd < maxBackwardRun)                
                    addSecondReversal(sm.segmentationClusters, st(j-1), st(j));
                    inreverse = false;
                else
                    addReversal(sm.segmentationClusters, st(j));
                    inreverse = true;
                end
        end
    end
    sc = sm.segmentationClusters;
    for j = 1:length(sc)
        sc(j).condenseKnownPoints;
    end
    kp = [sc(2:end).knownPoints];
    kp = kp([kp.track] == track);
    ri = track.isrun;
    ri([kp.inds]) = false;
    sc(1).addKnownPoints(track, find(ri)); %#ok<FNDSB>
    
    
end

                
function addOmegaTurn (sc, st, revst)
    
    existsAndDefault('revst', []);

    cr = st.track.getDerivedQuantity('scovRatio');
    dcr = st.track.getDerivedQuantity('dcovRatio');
    [mcr,I] = min(cr(st.inds));
    cI = st.centralInd;
    start = st.startInd;
    stop = st.endInd;
   
    if (start > length(cr) || stop > length(cr) || cI == length(cr))
        return;
    end

    curlstop = start - 1 + find(cr(start:(cI-1)) < 0.75*mcr + 0.25 * cr(start), 1, 'first');
    curlstart = cI + find(cr(cI+1):stop > 0.75*mcr + 0.25 * cr(stop), 1, 'first');
    if (isempty(curlstart) || isempty(curlstop) || curlstart > length(cr) || curlstop > length(cr))
        return;
    end
    stretchstart = max(curlstart+1, stop-2);
    curlupinds = start:curlstop;
    omegainds = (curlstop+1):(curlstart-1);
    uncurlinds = curlstart:(stretchstart - 1);
    stretchinds = stretchstart + (0:4);
    if (stretchinds(end) > length(cr))
        return;
    end
    names = {sc.name};
    sc(strcmpi('curl up', names)).addKnownPoints(st.track, curlupinds);
    sc(strcmpi('omega turn', names)).addKnownPoints(st.track, omegainds);
    sc(strcmpi('uncurl', names)).addKnownPoints(st.track, uncurlinds);
    sc(strcmpi('stretch', names)).addKnownPoints(st.track, stretchinds);
    
    if (~isempty(revst))
        sc(strcmpi('going backwards', names)).addKnownPoints(st.track, (revst.endInd+1):(start - 1));
    end
end

function addReversal(sc, st)
    dt = st.track.getDerivedQuantity('deltatheta');
    while (length(st.inds) < 5)
        st.startInd = max(1,st.startInd - 1);
        st.endInd = min(length(dt),st.endInd + 1);
        st.inds = st.startInd:st.endInd;
    end
    [~,I] = max(abs(dt(st.inds)));
    names = {sc.name};
    if (dt(st.startInd + I) > dt(st.startInd + I - 1))
        revinds = st.startInd + I + (-1:0);
    else
        revinds = st.startInd + I + (-2:-1);
    end
    intoinds = st.startInd:(revinds(1)-1);
    outofinds = (revinds(2)+1):st.endInd;
    sc(strcmpi('into reverse', names)).addKnownPoints(st.track, intoinds);
    sc(strcmpi('reverse', names)).addKnownPoints(st.track, revinds);
    sc(strcmpi('out of reverse', names)).addKnownPoints(st.track, outofinds);
end
   
function addSecondReversal (sc, prevst, st)
    dt = st.track.getDerivedQuantity('deltatheta');
    while (length(st.inds) < 5)
        st.startInd = max(1,st.startInd - 1);
        st.endInd = min(length(dt),st.endInd + 1);
        st.inds = st.startInd:st.endInd;
    end
    [~,I] = max(abs(dt(st.inds)));
    names = {sc.name};
    if (dt(st.startInd + I) > dt(st.startInd + I - 1))
        revinds = st.startInd + I + (-1:0);
    else
        revinds = st.startInd + I + (-2:-1);
    end
    intoinds = st.startInd:(revinds(1)-1);
    outofinds = (revinds(2)+1):st.endInd;
    sc(strcmpi('into second reversal', names)).addKnownPoints(st.track, intoinds);
    sc(strcmpi('second reversal', names)).addKnownPoints(st.track, revinds);
    sc(strcmpi('out of second reversal', names)).addKnownPoints(st.track, outofinds);
    
    sc(strcmpi('going backwards', names)).addKnownPoints(st.track, (prevst.endInd+1):(st.startInd - 1));
    
end
    
    
    
    