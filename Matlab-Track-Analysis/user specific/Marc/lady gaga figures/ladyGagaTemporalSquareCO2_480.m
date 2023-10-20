%run co2analysisTemporal20110110 first

existsAndDefault('loadCalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\square';
if (loadCalcs)
    load (fullfile (sd, 'co2_sq_480_calcs.mat'));%, 'co2_sq_ad', 'co2_sq_tc', 'co2_sq_tno');
    loadCalcs = false;
end

%}
labelnames = {'CO$_2$ 0-2.5\%', 'CO$_2$ 0-2.5\% gr63a'};
for j = 1:length(labelnames)
    po(j).legendEntry = labelnames{j};
    
    po(j).backgroundColor = 'w';
    po(j).onColor = 'k';
    po(j).alphaScale = 0.4;
    po(j).reverse = true;
end
sd = fullfile(sd, '8 min period');

whichGraphs = {'reorate_vs_time','speed_vs_time','reomag_vs_time','val_vs_time'};
existsAndDefault('savefiles', false);
if (savefiles)
    saveDirectory = fullfile(sd, 'with legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad_480,po,'showlegend', true, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs);

if (savefiles)
    saveDirectory = fullfile(sd, 'no legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad_480,po,'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', gcf + 1);
