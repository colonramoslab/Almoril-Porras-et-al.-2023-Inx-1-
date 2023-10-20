function [ respMat, tempMat, respCnt]  = respMaxFilter(fluor, temp, incFrames, respCnt)
%respFilter: Finds 'peaks' or responses and provides details of stim &
%response. Inputs are fluorescence and temperature time series with
%time-of-interest, in frames, as incFrames
%   respMat: Outputs 5sec window flanking responses. 
%   tempMat: Outputs 5sec stimulus flanking window.
%   respCnt:  Location of response within matrix

sampWin=50; % @10hz=5sec.

[ screen] = makeScreen( fluor.data, incFrames ); % make matrix mesh filter
[ respMat, tempMat, respCnt] = respMaxFind( fluor.data, temp.data,...
    'screen',screen,'sampWin',sampWin);
