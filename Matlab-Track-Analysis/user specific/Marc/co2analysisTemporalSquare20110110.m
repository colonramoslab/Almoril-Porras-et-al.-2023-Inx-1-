setCo2analysisFilesTemporalSquare

existsAndDefault('timerange', []);
existsAndDefault('trimrect', []);
%existsAndDefault('retrim', false(size(basedirs)));
existsAndDefault('resegment', false(size(basedirs)));
existsAndDefault('redogqs', false);

if (~exist('co2_eset_sq', 'var'))
 %   if (matlabpool('size') == 0)
  %      matlabpool
  %  end
    for j = 1:length(basedirs)
        co2_eset_sq(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
        if (isempty(co2_eset_sq(j).expt(end).globalQuantity))
            co2_eset_sq(j).executeExperimentFunction('addMetaDataFields');
           % co2_eset_sq(j).expt.addTonToff('
        end
        resegment(j) = true;
       
    end
%    matlabpool close
end
if (isempty(co2_eset_sq(end).expt(end).globalQuantity))
    for j = 1:length(co2_eset_sq)
        co2_eset_sq(j).executeExperimentFunction('addMetaDataFields');
        co2_eset_sq(j).expt.addTonToff('co2ppm1', 'square');
    end
end

for j = 1:length(co2_eset_sq)
%     if (retrim(j))
%         disp(['trimming - ' num2str(j)']);
%         co2_eset_sq(j).executeExperimentFunction('trimTracks', timerange, trimrect);
%         resegment(j) = true; %#ok<SAGROW>
%         retrim(j) = false;
%     end
     if resegment(j)
        disp(['segmenting - ' num2str(j)']);
%         t = [co2_eset_sq(j).expt.track];
%         for k = 1:length(t)
%             t(k).so.theta_cut = pi/2;
%         end
        co2_eset_sq(j).executeTrackFunction('setSegmentSpeeds');
        co2_eset_sq(j).executeTrackFunction('segmentTrack');
        resegment(j) = false; %#ok<SAGROW>
    end
    po(j).legendEntry = labelnames{j}; %#ok<SAGROW>
end
