function plotViterbi(sm, track, guessedStates, varargin)

c = 'rgcmyk'; c = [c c c c c];
msym = '******vvvvvv++++++xxxxxxhhhhhh';
if (isempty(guessedStates))
    guessedStates = sm.doViterbi(track);
end
valid = false(size(sm.segmentationClusters));
hold off;
for j = 1:length(sm.segmentationClusters)
    if (strcmpi(sm.segmentationClusters(j).name, 'run') || strcmpi(sm.segmentationClusters(j).name, 'runs'))
        marker = 'b.';
    else
        marker = [c(j) msym(j)];
    end
    if (any(guessedStates == j))
        track.plotPath('sloc', marker, 'inds', guessedStates == j, varargin{:});
        hold on;
        valid(j) = true;
    end
end
legend({sm.segmentationClusters(valid).name}, 'Location', 'Best');
