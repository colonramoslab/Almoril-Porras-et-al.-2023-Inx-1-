function sm = twelveStateWormModel(varargin)
%function sm = twelveStateWormModel(varargin)
%
%creates the structure for a segmentation model with runs, omega turns,
%and reversals, with transition states (no kinks)

sm = SegmentationModel.thirteenStateWormModel(varargin{:});
sm.segmentationClusters = sm.segmentationClusters([1,3:end]);
sm.allowedTransitions = sm.allowedTransitions([1,3:end],[1,3:end]);

 
 