function mraArr = fixExperiment(mra, expt, varargin)
%function mraArr = fixExperiment(mra, expt, varagin)
%detect segmentation problems, load track images, call visualize
%then reextract tracks, then unload loaded track images
%badlimit is the lower limit on the fraction of invalid HT detected
%default: 0.85;  pass 'badlimit', X to change
%passing 'reextract', false will abort resegmenting and will also
%not clear loaded track images

reextract = true;
badlimit = 0.85;
varargin = assignApplicable(varargin);
badlimit
existsAndDefault('mra', MaggotReAnalyzer);
%inds = expt.detectPossibleSegmentationProblems();
badratio = expt.evaluateTrackExpression('mean(track.getDerivedQuantity(''ihtValid''))');
inds = find(badratio <= badlimit);

mraArr = repmat(mra, size(inds));
clearim = false(size(inds));
for j = 1:length(inds)
    if (isempty(expt.track(inds(j)).pt(1).imData))
        try
            expt.reloadTrack(inds(j));
        catch me
            disp(me.getReport)
            expt
            inds(j)
        end
        clearim(j) = true;
    end
    mraArr(j) = mra.visualize(expt.track(inds(j)), 'hideAnalyze', true);
end

if (reextract)
    for j = 1:length(inds)
        mraArr(j).reExtractTrack([],varargin{:});
        if (clearim(j))
            [expt.track(inds(j)).pt.imData] = deal([]);
        end
    end
end