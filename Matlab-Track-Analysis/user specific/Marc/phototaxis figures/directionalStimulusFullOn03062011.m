existsAndDefault('loadeset', false);
if (loadeset)
    eset = ExperimentSet.fromMatFiles('E:\phototaxis\Extracted Phototaxis Data rdiff\Directional Gradients\All Ones\CS\matfiles\dir_45_cs');
    loadeset = false;
    eset.executeTrackFunction('setSegmentSpeeds');
    eset.executeTrackFunction('segmentTrack');
end

if (~exist('ad', 'var') && exist('eset', 'var'))
    sno = spatialNavigationMaggotAnalysis;
    sno.angleBinSize = 30;
    sno.reoBinSize = 30;
    sno.runTimeBinSize = 5;
    sno.minHS = 1;
    ad = spatialNavigationMaggotAnalysis(eset, sno);
    save('C:\Users\Marc\Documents\figures\liz poster phototaxis\directional experiment\full intensity\analyzedData_20110510.mat', 'ad', 'sno');
end


if (~exist('ad', 'var'))
    %load ('C:\Users\Marc\Documents\figures\liz poster phototaxis\directional experiment\full intensity\analyzedData.mat');
    load ('C:\Users\Marc\Documents\figures\liz poster phototaxis\directional experiment\full intensity\analyzedData_20110510.mat');

end

sd = 'C:\Users\Marc\Documents\figures\liz poster phototaxis\directional experiment\full intensity\';


tolightColor = [1 0.8 0.2];
todarkColor = [0.3 0.05 0.5];
upColor = [0.1 0.7 0.3];
downColor = [0 0.1 0.9];

po.directions = [-180,-90, 0, 90];
po.colors = {todarkColor, downColor, tolightColor, upColor};
po.legendEntry = {'away from light', 'south', 'into light', 'north'};

existsAndDefault('savefigs', false);
whichGraphs = {'ReoDirDistribution','RunLengthHistogram', 'DirectionHistogram', 'RunStartHistogram','ReorientationRateVsHeading', 'InstantaneousDeltaThetaVsTheta','SpeedVsDirection', 'ReoDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness', 'FirstHeadSwingHandedness'};
saveDirectory = [];
sf = 1;
if (savefigs)
    saveDirectory = fullfile(sd, 'no titles', 'with legend');
end
spatialNavigationMaggotFiguresSingleExpt(ad, sno, po, 'showtitle', false, 'SaveDirectory', saveDirectory, 'showlegend',true, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');

sf = sf + length(whichGraphs);
if (savefigs)
    saveDirectory = fullfile(sd, 'no titles', 'without legend');
end
spatialNavigationMaggotFiguresSingleExpt(ad, sno, po, 'showtitle', false, 'SaveDirectory', saveDirectory, 'showlegend',false,'startFigNum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');


sf = sf + length(whichGraphs);
if (savefigs)
    saveDirectory = fullfile(sd, 'with titles');
end
spatialNavigationMaggotFiguresSingleExpt(ad, sno, po, 'showtitle', true, 'showlegend',true,'SaveDirectory', saveDirectory, 'whichGraphs', whichGraphs,'startFigNum', sf, 'backgroundColor', 'w');
