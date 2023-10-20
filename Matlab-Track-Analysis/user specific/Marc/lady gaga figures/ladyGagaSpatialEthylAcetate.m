basedirs = {'E:\larvalco2\Extracted\Ethyl Acetate 233ppm CS\50mL air in 1.5 L air\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 2.25ppm CS\40mL air in 2 L air\CS', ...
    'E:\larvalco2\Extracted\Ethyl Acetate masking experiment\25ppm masking 25ppm in vlaves\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 331ppm Gr63a\50mL air in 1 L air\Gr63a',...
    'E:\larvalco2\Extracted\Ethyl Acetate pure Spatial Gr63a\50mL co2 in 2 L air\Gr63a', 'E:\larvalco2\Extracted\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS',...
    'E:\larvalco2\Extracted\Ethyl Acetate Or83b1\23ppm at output spatial\Or83b1', 'E:\larvalco2\Extracted\ethyl acetate Or83b2 spatial\23ppm middle in 2L air\Or83b2'}; 

esetnames = {'cs_ea_0_460', 'cs_ea_0_5', 'cs_ea_25_50', 'gr63a_ea_0_660', 'gr63a_ea_50_2', 'cs_ea_0_48', 'or83b_ea_0_46','Or83b2_ea_0_46'};

for j = 1:length(basedirs)
    d = dir(fullfile (basedirs{j}, 'matfiles', [esetnames{j} '_experiment*.mat']));
    if (isempty(d))
        error (['problem with ' num2str(j)']);
    end
end

%undiluted = 2200 ppm at 50 / 2000, which would say vapor pressure of ethyl
%acetate is 0.088 atm = 
%<math>\scriptstyle \log_{10} P_{mmHg} = 7.09808 - \frac {1238.71}
%{217.0+T}</math> = 82 mm mercury at 22 C = 0.108 atm -- close

%2000 x dilution at 50 / 2000 approx 20 ppm mean

descriptions = {'CS EA 20 ppm/cm', 'CS EA 0.2 ppm/cm', 'CS EA 2 ppm/cm, masked by 25 ppm EA', 'GR63a EA 30 ppm/cm', ...
    'GR63a EA 200 ppm/cm', 'CS EA 2 ppm/cm', 'OR83b1 EA 2 ppm/cm', 'OR83b2 EA 2 ppm/cm'};

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
     for j = 1:length(eset)
         eset(j).setAutocorrTau;
     end
     
end

if (~exist('ea_spatial_ad', 'var'))
    if (exist('eset', 'var'))
        ea_spatial_sno = spatialNavigationMaggotAnalysis;
        ea_spatial_sno.validname = 'eti';
        ea_spatial_sno.validoperation = @(x) x >= 120 & x <= 1020;
        ea_spatial_sno.minHS = 1;
        ea_spatial_sno.autocorr_timerange = [120 1020];
        ea_spatial_sno.runTimeBinSize = 5;
        ea_spatial_ad = spatialNavigationMaggotAnalysis(eset, ea_spatial_sno);
        if (savead)
            save('C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\ea_spatial_calcs.mat', 'ea_spatial_ad', 'ea_spatial_sno');
        end
        
    else
        load ('C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\ea_spatial_calcs.mat');
    end
end
ccc = 'kgcymbrm';
sss = 'sodvh>p^<';
for j = 1:length(descriptions)
    po(j).legendEntry = descriptions{j};
    po(j).color = ccc(j);
    po(j).marker = sss(j);
    
end
savedir = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\';
existsAndDefault('savefigs', false);

SaveDirectory = [];

if (savefigs)
    SaveDirectory = fullfile(savedir, 'all ea spatial');
     d = dir(SaveDirectory);
    if (isempty(d))
        mkdir(SaveDirectory);
    end
end

bc = 'w';
navorder = [1 6 2 3 7 8];
spatialNavigationMaggotFigures(ea_spatial_ad(navorder), ea_spatial_sno, po(navorder), 'whichGraphs', {'NavigationIndex'}, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);

if (savefigs)
    SaveDirectory = fullfile(savedir, 'cs ea 2 ppm per cm');
    d = dir(SaveDirectory);
    if (isempty(d))
        mkdir(SaveDirectory);
    end
end

startFigNum = 2;
toExamine = 6;
spatialNavigationMaggotFigures(ea_spatial_ad(toExamine), ea_spatial_sno, po(toExamine), 'whichGraphs', {'FirstHeadSwingBias','HeadSwingAcceptanceDirection'},'startFigNum', startFigNum, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);


%'DirectionHistogram',
%   'RunStartHistogram','ReorientationRateVsHeading',
%   'InstantaneousDeltaThetaVsTheta','SpeedVsDirection',
%   'ReoDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness',
%   'FirstHeadSwingHandedness', 'ReoDirDistribution'
startFigNum = 3;
whichGraphs = {'DirectionHistogram', 'ReorientationRateVsHeading',  'InstantaneousDeltaThetaVsTheta', 'SpeedVsDirection',  'ReoDirVsHeading','ReoMagVsHeading', 'ReoDirDistribution','HeadSwingAcceptanceHandedness','RunLengthHistogram','FirstHeadSwingHandedness'};
spatialNavigationMaggotFiguresSingleExpt(ea_spatial_ad(toExamine), ea_spatial_sno,[], 'whichGraphs', whichGraphs, 'startFigNum', startFigNum, 'backgroundColor', bc, 'SaveDirectory', SaveDirectory);
