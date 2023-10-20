basedirs = {'E:\larvalco2\Extracted\50 mL CO2 in 2 L air\CS', 'E:\larvalco2\Extracted\50 mL CO2 in 2 L air\gr63a',...
    'E:\larvalco2\Extracted\10mL CO2 2L air\CS', 'E:\larvalco2\Extracted\5p4 mL co2 in 2 L air\CS',...
    'E:\larvalco2\Extracted\control air 50 mL in 2L air\CS', 'E:\larvalco2\Extracted\no air flow\CS',...
    'E:\larvalco2\Extracted\20p4mL CO2 2L air\CS'};
esetnames = {'cs_50', 'gr63a_50', 'cs_10', 'cs_5', 'cs_0', 'cs_noair', 'cs_20'};
labelnames = {'CS 0-5\%', 'Gr63a 0-5\%', 'CS 0-1\%', 'CS 0-0.5\%', 'CS air control', 'CS no air control', 'CS 0-2\%'};
figdir = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 spatial';

for j = 1:length(basedirs)
    d = dir(fullfile (basedirs{j}, 'matfiles', [esetnames{j} '_experiment*.mat']));
    if (isempty(d))
        error (['problem with ' num2str(j)']);
    end
end

descriptions = {'CS CO$_2$ 2500 ppm/cm', 'Gr63a CO$_2$ 2500 ppm/cm', 'CS CO$_2$ 500 ppm/cm', 'CS CO$_2$ 250 ppm/cm', ...
    'CS clean air', 'CS no air flow', 'CS CO$_2$ 1000 ppm/cm'};

existsAndDefault('resegment', false);
existsAndDefault('loadeset', false);
existsAndDefault('savead', false);
if (loadeset && ~exist('eset', 'var'))
    ecl = ESetCleaner();
    ecl.rpmCut = 1;
    ecl.askFirst = false;
    ecl.showFigsInReport = false;
    for j = 1:length(basedirs)
        eset(j) = ExperimentSet.fromMatFiles(fullfile (basedirs{j}, 'matfiles', esetnames{j}));
        ecl.getReport(eset(j));
        ecl.clean(eset(j));
    end
    resegment = true;
end
if (exist('eset', 'var') && resegment)
     for j = 1:length(eset)
        eset(j).executeTrackFunction('setSegmentSpeeds');
        eset(j).executeTrackFunction('segmentTrack');
     end
     resegment = false;
     actual_tau = [eset.autocorr_tau];
     actual_tau
     [eset.autocorr_tau] = deal(25);
end

if (~exist('co2_spatial_ad', 'var'))
    if (exist('eset', 'var'))
        co2_spatial_sno = spatialNavigationMaggotAnalysis;
        co2_spatial_sno.validname = 'eti';
        co2_spatial_sno.validoperation = @(x) x >= 120 & x <= 1020;
        co2_spatial_sno.minHS = 1;
        co2_spatial_sno.autocorr_timerange = [120 1020];
        co2_spatial_sno.reoBinSize = 45;
        co2_spatial_sno.runTimeBinSize = 5;
        co2_spatial_ad = spatialNavigationMaggotAnalysis(eset, co2_spatial_sno);
        if (savead)
            save(fullfile(figdir, 'co2_spatial_calcs'), 'co2_spatial_ad', 'co2_spatial_sno');
        end
        
    else
        load (fullfile(figdir, 'co2_spatial_calcs'));
    end
end
ccc = 'kgcymbrm';
sss = 'sodvh>p^<';
for j = 1:length(descriptions)
    po(j).legendEntry = descriptions{j};
    po(j).color = ccc(j);
    po(j).marker = sss(j);
    
end

savedir = figdir;
existsAndDefault('savefigs', false);

SaveDirectory = [];

if (savefigs)
    SaveDirectory = fullfile(savedir, 'all co2 spatial');
    d = dir(SaveDirectory);
    if (isempty(d))
        mkdir(SaveDirectory);
    end
end

bc = 'w';

descriptions = {'CS CO$_2$ 250 ppm/cm', 'Gr63a CO$_2$ 250 ppm/cm', 'CS CO$_2$ 50 ppm/cm', 'CS CO$_2$ 25 ppm/cm', ...
    'CS clean air', 'CS no air flow', 'CS CO$_2$ 100 ppm/cm'};

navorder = [1 7 3 4 2 5 6];
spatialNavigationMaggotFigures(co2_spatial_ad(navorder), co2_spatial_sno, po(navorder), 'whichGraphs', {'NavigationIndex'}, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);

if (savefigs)
    SaveDirectory = fullfile(savedir, 'cs co2 2500 ppm per cm');
    d = dir(SaveDirectory);
    if (isempty(d))
        mkdir(SaveDirectory);
    end
end

startFigNum = 2;
toExamine = 1;
spatialNavigationMaggotFigures(co2_spatial_ad(toExamine), co2_spatial_sno, po(toExamine), 'whichGraphs', {'FirstHeadSwingBias','HeadSwingAcceptanceDirection'},'startFigNum', startFigNum, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);


%'DirectionHistogram',
%   'RunStartHistogram','ReorientationRateVsHeading',
%   'InstantaneousDeltaThetaVsTheta','SpeedVsDirection',
%   'ReoDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness',
%   'FirstHeadSwingHandedness', 'ReoDirDistribution'
startFigNum = gcf+1;
whichGraphs = {'DirectionHistogram', 'ReorientationRateVsHeading',  'InstantaneousDeltaThetaVsTheta', 'SpeedVsDirection',  'ReoDirVsHeading','ReoMagVsHeading', 'ReoDirDistribution','HeadSwingAcceptanceHandedness','RunLengthHistogram','FirstHeadSwingHandedness'};
spatialNavigationMaggotFiguresSingleExpt(co2_spatial_ad(toExamine), co2_spatial_sno,[], 'whichGraphs', whichGraphs, 'startFigNum', startFigNum, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);
