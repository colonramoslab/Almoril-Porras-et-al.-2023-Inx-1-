%run ethylAcetateTemporal20110124 first
%{
'Ethyl Acetate (conditioned) 0-1100 ppm, 220 ppm/min'
    'Ethyl Acetate (unconditioned) 0-2200 ppm, 440 ppm/min'
    'Ethyl Acetate (conditioned) 0-2200 ppm, 440 ppm/min'
    'Ethyl Acetate (conditioned) 1100-2200 ppm, 220 ppm/min'
    'Ethyl Acetate (unconditioned) 0-2200 ppm exponential'
    'Ethyl Acetate (unconditioned) 0-22 ppm exponential'
    'Ethyl Acetate (unconditioned) 0-0.2 ppm exponential'
    'Ethyl Acetate square 0-220 ppm'
    'Ethyl Acetate square 0-440 ppm'
    'Ethyl Acetate square 0-2200 ppm'
    'Ethyl Acetate Exponential 0 - 2 ppm'
%}
existsAndDefault('loadcalcs', true);
sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate temporal';
if (loadcalcs)
    load(fullfile(sd, 'ea_temp_calcs.mat'), 'tcsq', 'tctri', 'tcexp', 'tritno', 'exptno', 'sqtno', 'adtri', 'adexp', 'adsq');
    loadcalcs = false;
end
    
le = {'Ethyl Acetate 0-1100 ppm, 220 ppm/min',  'Ethyl Acetate 0-2200 ppm, 440 ppm/min (uc)', 'Ethyl Acetate 0-2200 ppm, 440 ppm/min', 'Ethyl Acetate 1100-2200 220 ppm/min',...
    'Ethyl Acetate exponential to 2200 ppm', 'Ethyl Acetate exponential to 22 ppm', 'Ethyl Acetate exponential to 0.2 ppm', 'Ethyl Acetate square 0-220', 'Ethyl Acetate square 0-440', 'Ethyl Acetate square 0-2200',...
    'Ethyl Acetate exponential to 2 ppm', 'Ethyl Acetate 500-2500 ppm, 400 ppm/min'};
triinds = [1:4 12];
expinds = [5:7 11];
sqinds = 8:10;

tritno = temporalNavigationAnalysis;
tritno.fieldname = 'vocppm';
tritno.period = 600;
tritno.minHS = 1;
tritno.timeBinSize = 30;
     
if (~exist ('tctri', 'var'))
    tctri = temporalCalculations(ea_eset(triinds), tritno);
end

exptno = temporalNavigationAnalysis;
exptno.fieldname = 'vocppm';
exptno.rampType = 'exponential';
exptno.period = 600;
exptno.minHS = 1;
exptno.timeBinSize = 30;
     
if (~exist ('tcexp', 'var'))
    tcexp = temporalCalculations(ea_eset(expinds), exptno);
end


sqtno = temporalNavigationAnalysis;
sqtno.fieldname = 'vocppm';
sqtno.rampType = 'square';
sqtno.period = 240;
sqtno.minHS = 1;
sqtno.timeBinSize = 12;
     
if (~exist ('tcsq', 'var'))
    tcsq = temporalCalculations(ea_eset(sqinds), sqtno);
end


if (~exist('adtri', 'var'))
    adtri = temporalNavigationAnalysis(ea_eset(triinds), tritno.fieldname, tritno, tctri);
end

if (~exist('adexp', 'var'))
    adexp = temporalNavigationAnalysis(ea_eset(expinds), exptno.fieldname, exptno, tcexp);
end
[adexp.rampType] = deal('exponential');

if (~exist('adsq', 'var'))
    adsq = temporalNavigationAnalysis(ea_eset(sqinds), sqtno.fieldname, sqtno, tcsq);
end
existsAndDefault('savecalcs', false);

sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate temporal';
if (savecalcs)
    save(fullfile(sd, 'ea_temp_calcs.mat'), 'tcsq', 'tctri', 'tcexp', 'tritno', 'exptno', 'sqtno', 'adtri', 'adexp', 'adsq');
    savecalcs = false;
end


existsAndDefault('savefiles', false);
for j = 1:length(le)
    po(j).legendEntry = le{j};
    po(j).backgroundColor = 'w';
    po(j).onColor = [0.6 0 0.6];
    po(j).alphaScale = 0.8;
    po(j).reverse = true;
end
whichGraphs = {'reorate_vs_time','speed_vs_time','reomag_vs_time','val_vs_time'};
if (savefiles)
    saveDirectory = fullfile(sd, 'triangle', 'with legends');
   
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
%whichtri = [1 3 4];
whichtri = 5;
temporalNavigationFigures(adtri(whichtri),po(triinds(whichtri)),'showlegend', true, 'SaveDirectory', saveDirectory, 'startFignum', 1, 'whichGraphs', whichGraphs);

if (savefiles)
    saveDirectory = fullfile(sd, 'exponential', 'with legends');
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
whichexp = [3 4 2 1];
temporalNavigationFigures(adexp(whichexp),po(expinds(whichexp)),'showlegend', true, 'SaveDirectory', saveDirectory, 'startFignum', gcf+1, 'whichGraphs', whichGraphs);


if (savefiles)
    saveDirectory = fullfile(sd, 'square', 'with legends');
    
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end
else
    saveDirectory = [];
end
temporalNavigationFigures(adsq,po(sqinds),'showlegend', true, 'SaveDirectory', saveDirectory, 'startFignum', gcf+1, 'whichGraphs', whichGraphs);

if (savefiles)
    saveDirectory = fullfile(sd, 'triangle', 'no legends');
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end

   % whichtri = [1 3 4];
    temporalNavigationFigures(adtri(whichtri),po(triinds(whichtri)),'showlegend', false, 'SaveDirectory', saveDirectory, 'startFignum', 1, 'whichGraphs', whichGraphs);
    
    saveDirectory = fullfile(sd, 'exponential', 'no legends');
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end

    whichexp = [3 4 2 1];
    temporalNavigationFigures(adexp(whichexp),po(expinds(whichexp)),'showlegend', false, 'SaveDirectory', saveDirectory, 'startFignum', gcf+1, 'whichGraphs', whichGraphs);
    
    
    saveDirectory = fullfile(sd, 'square', 'no legends');
    d = dir(saveDirectory);
    if (isempty(d))
        mkdir(saveDirectory);
    end

    temporalNavigationFigures(adsq,po(sqinds),'showlegend', false, 'SaveDirectory', saveDirectory, 'startFignum', gcf+1, 'whichGraphs', whichGraphs);
end
