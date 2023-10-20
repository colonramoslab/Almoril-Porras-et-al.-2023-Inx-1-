setCo2analysisFiles


if (~exist('co2_eset', 'var'))
    if (matlabpool('size') == 0)
        matlabpool
    end
    for j = 1:length(basedirs)
        co2_eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
    end
    matlabpool close
end

ecl = ESetCleaner;
ecl.minSpeed = 0.01;
ecl.minPts = 300;
ecl.showFigsInReport = false;
ecl.askFirst = false;
ts0 = tic;
for j = 1:length(co2_eset)
    ts = tic;
    disp ([num2str(j) '/' num2str(length(co2_eset)) ' - cleaning - ' num2str(toc(ts0)) ' s elapsed so far']); 
    ecl.clean(co2_eset(j));
    disp (['done. et = ' num2str(toc(ts)) ' calculating required quantities for segmentation']);
    dqs = {'curv', 'sspineTheta', 'vel_dp', 'speed'};
    for k = 1:length(dqs)
       co2_eset(j).gatherField(dqs{k});
    end
    disp (['done. et = ' num2str(toc(ts)) ' calculating segmentation speeds']);
    co2_eset(j).executeTrackFunction('setSegmentSpeeds');
    disp (['done. et = ' num2str(toc(ts)) ' saving eset to disk']);
    co2_eset(j).toMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j})); 
    disp (['done. total time = ' num2str(toc(ts))]);
end