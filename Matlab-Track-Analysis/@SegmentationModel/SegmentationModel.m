classdef SegmentationModel < handle
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        eset; % < ExperimentSet
        segmentationClusters; % < SegmentationCluster
        derivationRules;% rules by which parameters were derived
        weightOnKnown; %[0,1]; how much weight to place on the known points when doing EM
        allowedTransitions; %in a hmm, allowedTransitions(i,j) = true iff sc(i) can go to sc(j)
        hmm_prior; %prior probabilities, for hmm
        hmm_transmat; %transition probabilities for hmm
        hmm_mix; %mix ratio for hmm
        hmm_mean; %mean for clusters in hmm
        hmm_cov; %covariance matrix for clusters in hmm
        hmm_musub; %mean of all data, to subtract off for standardization
        hmm_sigmnorm; %sigma of all data, to subtract off for standardization
    end
    
    methods
        hm = toHashMap(sm,varargin);
        str = toString(sm,varargin);
        toFile(sm, filename, varargin);
        weightedEM(sm, varargin);
        bnet = createBNET (sm, varargin); %make a bnet hmm for this model
        cases = createBnetCases(sm, tracks, varargin); %create training data for this model
        initHMM (sm, varargin);
        [prior, transmat, mixmat, mu, sigma] = hmm_em(sm, tracks, varargin);
        [data, known_points] = createHMMData (sm, tracks, varargin);
        guessedStates = doViterbi (sm, track);
        plotViterbi(sm, track, guessedStates, varargin);
        sm2 = toSaveFriendly(sm, varargin);
        toMatFile (sm, filename);
        sm2 = fromSaveFriendly(sm, eset);
        setStandardization(sm, tracklist);
        addTransitionCluster (sm, basecluster, newclusterlocation, varargin);
        drawAllowedTransitions(sm);
        cnum = nextCluster(sm, currentcluster);
        cnum = prevCluster(sm, currentcluster);
        setDataFields(sm, fieldnames, operations);
    end
    
    methods(Static)
        sm = fromHashMap(hm, varargin);
        sm = fromString(str, varargin);
        sm = fromFile(fn, varargin);
        sm = fromMatFile(fn, eset);
        sm = fromSegmentedTracks(fn, src, varargin);
        sm = thirteenStateWormModel(varargin);
        sm = twelveStateWormModel(varargin);
    end
    
end

