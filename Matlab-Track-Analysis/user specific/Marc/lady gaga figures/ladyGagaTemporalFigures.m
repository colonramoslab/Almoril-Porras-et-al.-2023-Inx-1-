%run co2analysisTemporal20110110 first

existsAndDefault('loadCalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\figure 5\';
if (loadCalcs)
    load (fullfile ('C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\triangle', 'co2_tri_calcs.mat'), 'co2_tri_ad', 'co2_tri_tc', 'co2_tri_tno');
    load (fullfile ('C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\square', 'co2_sq_calcs.mat'), 'co2_sq_ad', 'co2_sq_tc', 'co2_sq_tno');
    load (fullfile ('C:\Users\Marc\Documents\figures\lady gaga paper\co2 temporal\square', 'co2_sq_480_calcs.mat'));
    loadCalcs = false;
end

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
clear po
po.backgroundColor = 'w';
po.onColor = [];
po.reverse = true;
po.shadedErrorRegion = true;
po = repmat(po, 10, 1);


whichtri = 1;
if (savefiles)
    saveDirectory = fullfile(sd, 'a');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_tri_ad(whichtri),po(whichtri),'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', 1);

if (savefiles)
    saveDirectory = fullfile(sd, 'b');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad, po,'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', gcf + 1);

if (savefiles)
    saveDirectory = fullfile(sd, 'c');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(co2_sq_ad_480, po,'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', gcf + 1);

