setCo2analysisFiles
ts10 = tic;
existsAndDefault('timerange', [120 1020]);
existsAndDefault('trimrect', []);
existsAndDefault('retrim', false(size(basedirs)));
existsAndDefault('resegment', false(size(basedirs)));

ecl = ESetCleaner;
ecl.askFirst = false;
ecl.minDist = 1.5;
ecl.minPts = 300;
if (~exist('co2_eset', 'var'))
    if (matlabpool('size') == 0)
        matlabpool
    end
    for j = 1:length(basedirs)
        co2_eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
        ecl.clean(co2_eset(j));
        resegment(j) = true;
    end
    matlabpool close
end
ccc = 'bgrcmyk';
sss = 'sodvh>p^<';

for j = 1:length(co2_eset)
    disp([num2str(j) '/' num2str(length(co2_eset))]);
    if (retrim(j))
        co2_eset(j).executeExperimentFunction('trimTracks', timerange, trimrect);
        resegment(j) = true; %#ok<SAGROW>
        retrim(j) = false;
    end
    if resegment(j)
        co2_eset(j).executeTrackFunction('setSegmentSpeeds');
        co2_eset(j).executeTrackFunction('segmentTrack');
        resegment(j) = false; %#ok<SAGROW>
    end
    po(j).legendEntry = labelnames{j}; %#ok<SAGROW>
    po(j).lineWidth = 2;
    po(j).color = ccc(mod(j-1, length(ccc)) + 1);
    po(j).marker = sss(mod(j-1, length(sss)) + 1);
    po(j).plotOptions = {};
end
sno.angleBinSize = 45;
sno.minHSTheta = 20;
sno.preferredDirection = 0;
if (~isempty(timerange))
    sno.validname = 'eti';
    sno.validoperation = @(x) logical(x >= timerange(1) & x <= timerange(2));
end
if (~exist ('ad', 'var'))
    ad = spatialNavigationMaggotAnalysis(co2_eset, sno);
end
varying_conc_inds = ([1 7 3 4 5]);
controls_inds = ([1 2 5 6]);

%{
ad_varyingconc = spatialNavigationMaggotAnalysis(co2_eset(varying_conc_inds), sno);
ad_controls = spatialNavigationMaggotAnalysis(co2_eset(controls_inds), sno);

spatialNavigationMaggotFigures(ad_varyingconc, sno, po(varying_conc_inds));
%}

spatialNavigationMaggotFigures(ad(varying_conc_inds), sno, po(varying_conc_inds));

toc(ts10)