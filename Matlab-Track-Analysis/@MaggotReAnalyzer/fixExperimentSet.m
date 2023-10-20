function mraArr = fixExperimentSet(mra, eset, varargin)
%function mraArr = fixExperimentSet(mra, eset, varargin)
%
%calls fix experiment on each experiment in the set
%then reextracts all together at end


existsAndDefault('mra', MaggotReAnalyzer);

for k = 1:length(eset.expt)
    mraArr{k} = mra.fixExperiment(eset.expt(k), 'reextract', false, varargin{:}); %#ok<AGROW>
    
end
disp('ok thank you:  I will now reextract tracks');
for k = 1:length(mraArr)
    for j = 1:length(mraArr{k})
        mraArr{k}(j).reExtractTrack([],varargin{:});
        [mraArr{k}(j).track.pt.imData] = deal([]); %#ok<AGROW>
    end
end