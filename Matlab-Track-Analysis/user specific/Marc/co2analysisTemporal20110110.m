setCo2analysisFilesTemporal

existsAndDefault('timerange', []);
%existsAndDefault('trimrect', []);
%existsAndDefault('retrim', false(size(basedirs)));
existsAndDefault('resegment', false(size(basedirs)));
existsAndDefault('redogqs', true);
if (~exist('co2_eset', 'var'))
    if (matlabpool('size') == 0)
        matlabpool
    end
    for j = 1:length(basedirs)
        co2_eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
        resegment(j) = true;
        %retrim(j) = false;
        if (redogqs || isempty(co2_eset(j).expt(end).globalQuantity))
            co2_eset(j).executeExperimentFunction('addMetaDataFields');
            co2_eset(j).expt.addTonToff('co2ppm1', 'triangle');
        end
    end
    matlabpool close
end
%{
if (isempty(co2_eset(end).expt(end).globalQuantity))
    for j = 1:length(co2_eset)
        co2_eset(j).executeExperimentFunction('addMetaDataFields');
        co2_eset(j).expt.addTonToff('co2ppm1', 'triangle');
    end
end
%}
for j = 1:length(co2_eset)
    %{
    if (retrim(j))
        disp(['trimming - ' num2str(j)']);
        co2_eset(j).executeExperimentFunction('trimTracks', timerange, trimrect);
        resegment(j) = true; %#ok<SAGROW>
        retrim(j) = false;
    end
    %}
    if resegment(j)
        disp(['segmenting - ' num2str(j)']);
        co2_eset(j).executeTrackFunction('setSegmentSpeeds');
        co2_eset(j).executeTrackFunction('segmentTrack');
        resegment(j) = false; %#ok<SAGROW>
    end
    po(j).legendEntry = labelnames{j}; %#ok<SAGROW>
end
