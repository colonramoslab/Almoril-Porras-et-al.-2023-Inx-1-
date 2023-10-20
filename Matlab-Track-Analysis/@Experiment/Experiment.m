%Experiment
%A set of tracks extracted from the same experiment
%Each experiment is loaded from a single .bin file created by the
%track extraction software


classdef Experiment < handle
    %a set of tracks extracted from the same experiment
    
    properties
        fname = ''; %the name of the .bin file representing the experiment
        timfname = ''; %the name of the timing file that tells when each frame happened
        camcalinfo = []; %camera calibration information; not tested at the moment
        elapsedTime = []; %the elapsed time from the start of the experiment for each frame
        temperature = []; %???
        
        dr = DerivationRules(); %rules for interpolating, smoothing, and differentiating
        so = WormSegmentOptions(); %rules for segmenting tracks into runs & reorientations
        globalQuantity; %used to add extra information (e.g. light, temperature) to the tracks
    end
    properties (AbortSet = true)
        track; %the set of tracks contained in the experiment; see Track, MaggotTrack
    end
    
    properties (Transient = true)
        fid = 0;
    end
    
    methods
        addtime (expt, timfname)
        addGlobalQuantity(expt,varargin);
        assignGlobalQuantities(expt,varargin);
        addStandardizedField(expt, fieldname, varargin);
        addMetaDataFields(expt, timfname, varargin);
        infertime (expt, fstub)
        stitchTracks (expt, maxFrameInterval, maxDist, varargin)
        cleanTracks (expt, minFrames, minDist)
        qvec = gatherField(expt, fieldname, varargin)
        qvec = gatherFieldInTracks(expt, fieldname, tracknums, varargin)
        qvec = gatherFromSubField(expt, subfield, fieldname, varargin);
        qvec = gatherSubField (expt, field, subfield, varargin)
        [pt, track, trackind, ptind] = findNearestPoint (expt, x, y);
        pt = reloadPoint (expt, pt);
        track2 = reloadTrack(expt, track);
        im = pointIm (expt, varargin);
        openDataFile(expt);
        closeDataFile(expt);
        calculateDerivedQuantity(expt, quantityNames, recalculate);
        flagOmegaTurnsAndReversals(expt);
        flagReorientations(expt);
        segmentTracks(expt, segmentOptions);
        varargout = executeTrackFunction(expt, func, varargin);
        result = evaluateTrackExpression(expt, expression);
        detectCollisions(expt, maxdist);
        im = makeWholeFrameImage (expt, ind, varargin);
        
        %trimTracks cuts out any part of the track outside
        %min(timerange),max(timerange) and removes any part of the track
        %from the point the track leaves validrect until the end of the
        %track; leave timerange or validrect empty to disable
        trimTracks(expt, timerange, validrect);
        
        %pruneTracks removes completely any track that starts outside
        %starttimerange in elapsedTime, as well as any track that starts
        %outside startrect; leave startttimerange or startrect empty to disable
        pruneTracks(expt, starttimerange, startrect);
        
        [ps, f] = powerSpectrum(expt, quantityName, timeInterval, varargin);
        [xc, np, tx, nt] = crosscorrelate (expt, fieldname1, fieldname2, varargin);
        [ac, np, tx, nt] = autocorrelate (expt, fieldname, varargin);
        updateTrackSegmentationOptions(expt);
        setTimingByFrameDifference(expt, deltaT, overrideExisting);
        
        addTonToff(expt, fieldname, ramptype, varargin);
    end
    
    methods %set methods
        function set.track(obj, value)
            if (~isempty(value) && isa (value, 'Track'))
                [value.expt] = deal(obj);
            end
            obj.track = value;
        end
        
    end
    
    methods %constructor
        function expt = Experiment(varargin)
            if (~isempty(varargin) && isa(varargin{1}, 'Experiment'))
                fn = intersect(fieldnames(varargin{1}), fieldnames(expt));
                for j = 1:length(fn)
                    expt.(fn{j}) = varargin{1}.(fn{j});
                end
            end
        end
    end
    
    
    methods(Static)
        expt = fromFile (fname, timfname, loadContour, camcalinfo, minTrackLength)
    end
end

