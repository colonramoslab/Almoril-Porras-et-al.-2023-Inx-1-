function mra = visualize(mra,track,varargin)
%function mra = visualize(mra,track,varargin)
%
%assigns track to mra (optional)
%then calls the GUI to manipulate mra
%checks to make sure that there is a track before proceeding

existsAndDefault('track', mra.track);
mra.track = track;
if(isempty(mra.track))
    disp ('you must provide a track in call or in mra to proceed');
    return;
end
if (isempty(mra.targetArea))
    mra.targetArea = track.pt(1).targetArea;
end
mra = TrackReanalyzerGui('mra', mra,varargin{:});