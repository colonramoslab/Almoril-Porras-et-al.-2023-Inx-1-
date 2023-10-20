%run co2analysisTemporal20110110 first

existsAndDefault('loadCalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\square\revised';
if (loadCalcs)
    load (fullfile (sd, 'co2_sq_calcs.mat'), 'co2_sq_ad', 'co2_sq_tc', 'co2_sq_tno');
    loadCalcs = false;
end
if (~exist ('co2_sq_tno', 'var'))
    co2_sq_tno = temporalNavigationAnalysis;
    co2_sq_tno.rampType = 'square';
    co2_sq_tno.fieldname = 'co2ppm1';
    co2_sq_tno.period = 240;
    co2_sq_tno.minHS = 0;
    co2_sq_tno.timeBinSize = 24;%changed from 20 to 30 to see effect -- was 10 for figures pre 5/26

    co2_sq_tno = repmat(co2_sq_tno, length(co2_eset_sq), 1);
end

existsAndDefault('savefiles', false);
if (~exist('co2_sq_tc', 'var'))
    for j = 1:3
        co2_sq_tc(j) = temporalCalculations(co2_eset_sq(j), co2_sq_tno(j));
    end
end

if (~exist('co2_sq_ad', 'var'))
    for j = 1:3
        co2_sq_ad(j) = temporalNavigationAnalysis(co2_eset_sq(j), co2_sq_tno(j).fieldname, co2_sq_tno(j), co2_sq_tc(j));
    end
end

%{
if (~exist('adlgbins', 'var'))
    for j = 1:4
        adlgbins(j) = temporalNavigationAnalysis(co2_eset_sq(j), co2_sq_tnolgbins(j).fieldname, co2_sq_tnolgbins(j), tc(j));
    end
end
%}
labelnames = {'CO$_2$ 0-0.25\%', 'CO$_2$ 0-.5\%', 'CO$_2$ 0-2.5\%'};
for j = 1:length(labelnames)
    po(j).legendEntry = labelnames{j};
    
    po(j).backgroundColor = 'w';
    po(j).onColor = [];
    po(j).alphaScale = 0.4;
    po(j).reverse = true;
    po(j).shadedErrorRegion = true;
end

existsAndDefault('saveCalcs', false);
if (saveCalcs)
    save (fullfile (sd, 'co2_sq_calcs.mat'), 'co2_sq_ad', 'co2_sq_tc', 'co2_sq_tno');
    saveCalcs = false;
end

whichGraphs = {'reorate_vs_time','speed_vs_time','reomag_vs_time','val_vs_time'};
%{
if (savefiles)
    saveDirectory = fullfile(sd, 'with legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad,po,'showlegend', true, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs);
%}
if (savefiles)
    saveDirectory = fullfile(sd, 'no legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad,po,'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', 1);
