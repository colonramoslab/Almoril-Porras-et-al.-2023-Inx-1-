classdef SegmentationCluster < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name; %descriptive name of this cluster
        isSharpTurn; %whether this cluster represents a sharp turn
        typeCode; %the type code of the sharp turn if isSharpTurn
        priorProbability; %the total fraction of time the animal spends 
                          %doing this behavior
        datafields; %which fields are used to cluster; K fields
        operation; %operations (e.g. abs, square, sin) on the fields; K operations
        clustMean; %cluster mean for each datafield; Kx1
        clustCov; %covariance for all datafields; KxK
        knownPoints; %Npts vector of kp struct;  kp has fields track, ind
                     %where the ind(s) (interpolated) pt(s) of kp.track
                     %is/are known to be a member of this cluster
    end
    
    methods
        lk = likelihood (sc, track);
        postprob = posteriorProbabilities(scSet, track);
        addKnownPoints(sc, track, inds); 
        calculateBasedOnKnownPoints(sc);
        [cm, ccv, cmr] = meanAndCovOfKnownPoints(sc, numclusters, varargin);
        [cm, cov] = meanAndCovEMStep(sc, point_source, postprob, varargin);
        str = toString(sc);
        toFile(sc, fname);
        hm = toHashMap(sc, varargin);
        [data, mu, sigma2] = getData(sc, point_source, varargin);
        sc2 = toSaveFriendly(sc);
        sc2 = fromSaveFriendly(sc, eset);
        toMatFile(sc, fname);
        condenseKnownPoints(sc);
        
    end
    methods (Static)
        sc = fromHashMap(obj, varargin);
        sc = fromString(str, varargin);
        sc = fromFile(fname, varargin);
        sc = fromMatFile (fname, eset);
    end
end

