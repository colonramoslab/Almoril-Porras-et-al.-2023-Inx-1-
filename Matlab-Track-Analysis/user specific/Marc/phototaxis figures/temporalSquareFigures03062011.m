if (~exist('tempad_5', 'var'))
    load ('E:\phototaxis\Extracted Phototaxis Data rdiff\Temporal Gradients\Square\CantonS\200s\2nds\temporalanalysis.mat');
end

sd = ('C:\Users\Marc\Documents\figures\liz poster phototaxis\temporal square wave');

existsAndDefault('savefigs', false);

po.offColor = [0.2 0.2 0.2];
%po.offColor = [0.1 0.1 0.1];
po.onColor = [];

whichGraphs = {'reorate_vs_time','reomag_vs_time','numhs_vs_time'};
saveDirectory = [];
sf = 1;
if (savefigs)
    saveDirectory = fullfile(sd, 'no titles');
end
temporalNavigationFigures(tempad_5, po, 'showtitle', false, 'SaveDirectory', saveDirectory, 'whichGraphs', whichGraphs);

sf = sf + length(whichGraphs);
if (savefigs)
    saveDirectory = fullfile(sd, 'with titles');
end
temporalNavigationFigures(tempad_5, po, 'showtitle', true, 'SaveDirectory', saveDirectory, 'whichGraphs', whichGraphs,'startFignum', sf);

po.offColor = [];
sf = sf + length(whichGraphs);
if (savefigs)
    saveDirectory = fullfile(sd, 'no backgrounds no titles');
    mkdir(saveDirectory);
end
temporalNavigationFigures(tempad_5, po, 'showtitle', false, 'SaveDirectory', saveDirectory, 'whichGraphs', whichGraphs,'startFignum', sf);

sf = sf + length(whichGraphs);
if (savefigs)
    saveDirectory = fullfile(sd, 'no backgrounds with titles');
    mkdir(saveDirectory);
end
temporalNavigationFigures(tempad_5, po, 'showtitle', true, 'SaveDirectory', saveDirectory, 'whichGraphs', whichGraphs,'startFignum', sf);
