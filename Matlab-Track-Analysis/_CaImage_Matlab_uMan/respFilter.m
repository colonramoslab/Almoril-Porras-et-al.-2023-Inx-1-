function [ fluorMat, tempMat, respCnt, peakMat]  = respFilter(fluor, temp, incFrames, varargin)
%respFilter: Finds 'peaks' or responses and provides details of stim &
%response. Inputs are fluorescence and temperature time series with
%time-of-interest, in frames, as incFrames
%   respMat: Outputs 5sec window flanking responses. 
%   tempMat: Outputs 5sec stimulus flanking window.
%   respCnt:  Location of response within matrix

sampWin=50; % @10hz=5sec.
derThresh=0.01; % 0.015 for AIY (190604), 0.002 for AFD (190604

%% get any input values
varargin = assignApplicable (varargin);


%% do some things...
[ screen] = makeScreen( fluor.data, incFrames ); % make matrix mesh filter

% [ fluorMat, tempMat, respCnt, peakMat] = respFind( fluor.data,
% temp.data,... % SWITCH from respFind to new respFind2 on 200420

[ fluorMat, tempMat, respCnt, peakMat] = respFind2( fluor.data, temp.data,...
    'screen',screen,'sampWin',sampWin, 'derThresh', derThresh);

end

