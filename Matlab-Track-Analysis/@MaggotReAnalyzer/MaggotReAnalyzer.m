classdef MaggotReAnalyzer
    %Re-analyzes an already loaded track of either ImTrackPoints or
    %TrackPoints
    %Purpose is to allow flexibility and correct errors from high speed c++
    %analysis program
    
    properties
        debug = false;
        track;
        targetArea;
        contourScale = 1;
        maxContourAngle = deg2rad(90);
    end
    
    methods
        %newpt = rethreshold (mra, pt, varargin)
        %rethresholds the current point to produce
        %a new MaggotTrackPoint
        %updates the contour threshold targetArea area loc
        newpt = rethreshold (mra, pt, varargin);
        
        %pt = findHT (mra, pt, varargin);
        %given a point with contour computed
        %finds head & tail candiates and saves one as pt.head
        %and the other as pt.tail
        %does not align head and tail
        pt = findHT (mra, pt, varargin);
        
        %findHT (mra, track, pt, varargin)
        %pt = alignHTPt (mra, pt, prevpt, varargin)
        %aligns the head and tail position of point 
        %with prevPt as long as both have valid HT
        
        pt = alignHTPt (mra, pt, prevpt, varargin);
        
        %alignHTSegment (mra, track, inds)
        %aligns the head/tail direction to the direction of motion
        %for a segment
        alignHTSegment (mra, track, inds)
        
        %alignHTTrack (mra, track)
        %aligns the head/tail direction to the direction of motion
        %for all segments of valid HT locations in a track
        alignHTTrack (mra, track)
        
        %reExtractTrack(mra, track, varargin)
        %
        %rethresholds, finds the head and tail, and aligns to the direction
        %of motion
        reExtractTrack(mra, track, varargin)
        
        %redoHTTrack(mra, track, varargin)
        %
        %from contour as computed, finds a refined head and tail location
        redoHTTrack(mra, track, varargin)
        
        
        %visualize(mra)
        %
        %visually set reanalysis params using GUI
        mra = visualize(mra,track,varargin);
        
        %fixExperiment
        %
        %detect segmentation problems, reload tracks, call visualize
        %then resegment all
        mraArr = fixExperiment(mra, expt, varargin);
        
        %calls fix experiment on each experiment in the set
        %then reextracts all together at end
        mraArr = fixExperimentSet(mra, eset, varargin);
        
        
    end
    
    methods (Static)
        %GUI
    end
    
end

