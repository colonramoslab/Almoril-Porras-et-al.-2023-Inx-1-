%all bin files will be sorted by date,
% clear cond;
% 
% j = 1;
% cond(j).name = 'low intensity directional proj high cs';
% cond(j).basedir = 'e:\Phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Ones\CS';
% cond(j).calibration = {'gershow2010-11-15a', 'Kane2010-12-20'};
% cond(j).whichcalibration = [1 1 1 1 2 2 2 2 2];
% cond(j).esetname = 'csDir45mid';
% 
% j = j+1;
% cond(j).name = 'low intensity directional proj high rh5rh6';
% cond(j).basedir = 'e:\Phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Ones\Rh5Rh6';
% cond(j).calibration = {'gershow2010-11-15a', 'Kane2010-12-20'};
% cond(j).whichcalibration = [1 2 2 2 2 2];
% cond(j).esetname = 'rh5rh6Dir45mid';
% 
% j = j+1;
% cond(j).name = 'low intensity directional proj low cs';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Zeros';
% cond(j).calibration = {'Kane2010-12-20', 'Kane2011-02-09'};
% cond(j).whichcalibration = [1 1 2 2 2 2];
% cond(j).esetname = 'csDir45low';
% 
% j = j+1;
% cond(j).name = 'high intensity directional cs';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\CantonS';
% cond(j).calibration = {'Liz2011-03-31'};
% cond(j).whichcalibration = [1 1 1 1];
% cond(j).esetname = 'wcsDir45high';
% 
% 
% j = j+1;
% cond(j).name = 'high intensity directional wcs';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\wCs';
% cond(j).calibration = {'Liz2011-03-31'};
% cond(j).whichcalibration = [1 1 1];
% cond(j).esetname = 'wcsDir45high';
% 
% 
% j = j+1;
% cond(j).name = 'high intensity directional gmr-hid';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\GMR-hid';
% cond(j).calibration = {'Liz2011-03-31'};
% cond(j).whichcalibration = [1 1 1 1];
% cond(j).esetname = 'gmrhidDir45high';
% 
% j = j+1;
% cond(j).name = 'high intensity directional rh5rh6';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\ywrh5rh6';
% cond(j).calibration = {'Liz2011-03-31'};
% cond(j).whichcalibration = [1 1 1 1];
% cond(j).esetname = 'rh5rh6Dir45high';
% 
% j = j+1;
% cond(j).name = 'projector off control';
% cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\Projector Off Control';
% cond(j).calibration = {'Kane2011-02-02'};
% cond(j).whichcalibration = [1 1 1 1];
% cond(j).esetname = 'csDirProfOff';
ts11 = tic;
drevli_setdirectionalfiles
%process camera calibrations
caldir = 'E:\phototaxis\calibrations';
existsAndDefault('redocameracalibrations', false);
disp ('processing camera calibrations');

for j = 1:length(cond)
    for k = 1:length(cond(j).calibration)
        cname = cond(j).calibration{k};
        cname(cname=='-') = '_';
        cname = ['cc_' cname]; %#ok<AGROW>
        d = dir(fullfile(cond(j).basedir, [cname '.mat']));
        if (~isempty(d) && ~redocameracalibrations)
            continue;
        end
        im = imread(fullfile(caldir, cond(j).calibration{k}, 'IR checker.tiff'));
        cc = CameraCalibration(im, 'flatten', true, 'flipy', true);
        eval([cname ' = cc;']);
        save(fullfile(cond(j).basedir, [cname '.mat']), cname);
    end
    close all
    redocameracalibrations = false;
end
            
disp ('done with camera calibrations');
toc(ts11)

%process bins to matfiles
timerange = [300 Inf];
ecl = ESetCleaner;
ecl.minHTValid = 0.9;
ecl.minSpeed = 0.5/60;
ecl.minDist = 1;
ecl.rpmCut = 2;
ecl.minPts = 600;

existsAndDefault('processbins', false);
if (processbins)
    disp('processing bins');
    for j = 1:length(cond)
        %set up the correct camera calibrations
        clear cname;
        for k = 1:length(cond(j).calibration)
            cname{k} = cond(j).calibration{k};
            cname{k}(cname{k}=='-') = '_';
            cname{k} = ['cc_' cname{k}];
            load (fullfile(cond(j).basedir, [cname{k} '.mat']), cname{k});
        end
        clear cc;
        for k = 1:length(cond(j).whichcalibration)
            eval(['cc(k) = ' cname{cond(j).whichcalibration(k)} ';']);
        end
        
        eset(j) = ExperimentSet.loadTrimStitchAndSave(cond(j).basedir, cond(j).esetname, ecl, cc, 'sortbydate', true); %#ok<SAGROW>
    end
    disp ('done processing bins');
    toc(ts11)
    processbins = false;
    disp ('segmenting tracks');
    for j = 1:length(eset)
        eset(j).executeTrackFunction('segmentTrack');
    end
    disp ('segmented!');
    toc(ts11);
end

%if we didn't process the bins, load the processed .matfiles
if (~exist('eset', 'var'))
    disp ('loading eset');
    for j = 1:length(cond)
        disp(['loading ' num2str(j) ' / ' num2str(length(cond))]);
        eset(j) = ExperimentSet.fromMatFiles(fullfile(cond(j).basedir, 'matfiles', cond(j).esetname)); %#ok<SAGROW>
        eset(j).executeTrackFunction('segmentTrack');
    end
    disp ('done loading eset');
    toc(ts11);
end
%{ 
%current defaults

          angleBinSize: 30
        dtAngleBinSize: 20
      hsdtAngleBinSize: 20
    preferredDirection: 0
            reoBinSize: 30
             hsBinSize: 90
          hsBinSpacing: 90
                 minHS: 0
            minHSTheta: 20
             validname: []
      relativeDirField: []
        dirOffsetField: []
        validoperation: @(x)logical(setNonFiniteToZero(x))
       confidenceLevel: 0.9500
    autocorr_timerange: []
        runTimeBinSize: 10
%}
calcdir = 'e:\phototaxis\dr evli\calculations\directional';
d = dir(calcdir);
if (isempty(d))
    mkdir(calcdir);
end

sno.minHS = 1;
disp('analyzing data')
ts12 = tic;
ad_dir = spatialNavigationMaggotAnalysis(eset, sno);
toc(ts12)
% disp ('testing save')
% add = ad_dir(1);
% fn = fieldnames(add);
% for j = 1:length(fn)
%     disp(fn{j});
%     ts99 = tic;
%     adtest.(fn{j}) = add.(fn{j});
%     save(fullfile(calcdir, 'test_save_directional_calculations.mat'), 'adtest');
%     toc(ts99)
% end
    

save(fullfile(calcdir, 'directional_calculations.mat'), 'ad_dir');
disp (['calculations took ' num2str(toc(ts12))]);
toc(ts11)
