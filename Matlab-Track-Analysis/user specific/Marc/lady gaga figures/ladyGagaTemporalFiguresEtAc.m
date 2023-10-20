%run co2analysisTemporal20110110 first

existsAndDefault('loadcalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\figure 6\';

if (loadcalcs)
    load(fullfile( 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate temporal', 'ea_temp_calcs.mat'), 'tcsq', 'tctri', 'tcexp', 'tritno', 'exptno', 'sqtno', 'adtri', 'adexp', 'adsq');
    loadcalcs = false;
end

whichGraphs = {'reorate_vs_time','speed_vs_time','reomag_vs_time','val_vs_time'};

existsAndDefault('savefiles', false);

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
whichtri = 5;
temporalNavigationFigures(adtri(whichtri),po(whichtri),'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', 1);

if (savefiles)
    saveDirectory = fullfile(sd, 'b');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(adsq, po,'showlegend', false, 'SaveDirectory',saveDirectory, 'whichGraphs', whichGraphs, 'startFignum', gcf + 1);

