%run co2analysisTemporal20110110 first

existsAndDefault('loadCalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\triangle';
if (loadCalcs)
    load (fullfile (sd, 'co2_tri_calcs.mat'), 'co2_tri_ad', 'co2_tri_tc', 'co2_tri_tno');
    loadCalcs = false;
end

if (~exist('co2_tri_tno', 'var'))
    co2_tri_tno = temporalNavigationAnalysis;
    co2_tri_tno.fieldname = 'co2ppm1';
    co2_tri_tno.period = 600;
    co2_tri_tno.minHS = 1;
    co2_tri_tno.timeBinSize = 30;

    co2_tri_tno = repmat(co2_tri_tno, 4, 1);
    co2_tri_tno(4).fieldname = 'co2ppm0';
end
existsAndDefault('savefiles', false);


if (~exist('co2_tri_tc', 'var'))
    for j = 1:4
        co2_tri_tc(j) = temporalCalculations(co2_eset(j), co2_tri_tno(j));
    end
end

if (~exist('co2_tri_ad', 'var'))
    for j = 1:4
        co2_tri_ad(j) = temporalNavigationAnalysis(co2_eset(j), co2_tri_tno(j).fieldname, co2_tri_tno(j), co2_tri_tc(j));
    end
end

%{
if (~exist('adlgbins', 'var'))
    for j = 1:4
        adlgbins(j) = temporalNavigationAnalysis(co2_eset(j), co2_tri_tnolgbins(j).fieldname, co2_tri_tnolgbins(j), tc(j));
    end
end
%}
le = {'CO$_2$ 0-5\%, 1\%/min', 'CO$_2$ 0-2.5\%, 0.5\%/min', 'CO$_2$ 2.5-5\%, 0.5\%/min','gr63a CO$_2$ 0-5\%, 1\%/min'};

for j = 1:4
    po(j).legendEntry = le{j};
    
    po(j).backgroundColor = 'w';
    %po(j).onColor = 'k';
    po(j).onColor = [];
    po(j).shadedErrorRegion = true;
    po(j).alphaScale = 0.4;
    po(j).reverse = true;
end

existsAndDefault('saveCalcs', false);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\triangle';
if (saveCalcs)
    save (fullfile (sd, 'co2_tri_calcs.mat'), 'co2_tri_ad', 'co2_tri_tc', 'co2_tri_tno');
    saveCalcs = false;
end
whichGraphs = {'reorate_vs_time','speed_vs_time','reomag_vs_time','val_vs_time'};

if (savefiles)
    saveDirectory = fullfile(sd, 'with legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
whichtri = [2 1 3];
temporalNavigationFigures(co2_tri_ad(whichtri),po(whichtri),'showlegend', true, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs);

if (savefiles)
    saveDirectory = fullfile(sd, 'no legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_tri_ad(whichtri),po(whichtri),'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', gcf + 1);

