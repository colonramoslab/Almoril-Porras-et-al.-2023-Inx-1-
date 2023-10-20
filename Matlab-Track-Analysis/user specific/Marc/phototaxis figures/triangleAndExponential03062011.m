if (~exist('ad_tri_20', 'var'))
    load ('C:\Users\Marc\Documents\figures\liz poster phototaxis\temporal exponential and triangle\triangle\analyzedData400.mat');
end
if (~exist('exp_ad_20', 'var'))
    load ('C:\Users\Marc\Documents\figures\liz poster phototaxis\temporal exponential and triangle\exponential\analyzedData400.mat');
end
sd = 'C:\Users\Marc\Documents\figures\liz poster phototaxis\temporal exponential and triangle';
existsAndDefault('savefigs', false);

lumFromCamExp = @(x) 30*x;
lumFromCamLin = @(x) 46*x;

llet = -200:200;
lltri = lumFromCamLin(200-abs(llet) + 1);

llexp = lumFromCamExp(.9961 * exp(0.025*(200-abs(llet))));

dlltri = 0.5*(diff(lltri([end 1:end])) + diff(lltri([1:end 1])));
dllexp = 0.5*(diff(llexp([end 1:end])) + diff(llexp([1:end 1])));

figure(1); clf(1); 
set(gcf, 'Color', 'w');
plot (llet, lltri, 'b', llet, llexp, 'g', 'LineWidth', 2); 
title ('Light Level vs. Time in Cycle');
ylabel ('light level (lux)');
xlabel ('time in cycle');
legend ('triangle', 'exponential', 'Location', 'BestOutside');


figure(2); clf(2); 
set(gcf, 'Color', 'w');
plot (llet, lltri/max(lltri), 'b', llet, llexp/max(llexp), 'g', 'LineWidth', 2); 
title ('Light Level vs. Time in Cycle');
ylabel ('light level / max light level');
xlabel ('time in cycle');
legend ('triangle', 'exponential', 'Location', 'BestOutside');

figure(3); clf(3);
set(gcf, 'Color', 'w');
plot (llet, dlltri, 'b', llet, dllexp, 'g', 'LineWidth', 2);
title ('Light Level Change vs. Time in Cycle');
ylabel ('$\Delta$ light level (lux/s)');
xlabel ('time in cycle');
legend ('triangle', 'exponential', 'Location', 'BestOutside');


figure(4); clf(4);
set(gcf, 'Color', 'w');
plot (llet, dlltri./lltri, 'b', llet, dllexp./llexp, 'g', 'LineWidth', 2); ylim([-0.1 0.1])
title ('Fractional Light Level Change vs. Time in Cycle');
ylabel ('$\Delta$ light level/light level');
xlabel ('time in cycle');
legend ('triangle', 'exponential', 'Location', 'BestOutside');
SaveDirectory = [];
if (savefigs)
    SaveDirectory = fullfile(sd, 'combined');
end
saveName = {'light vs time absolute', 'light vs time scaled', 'change vs time', 'fractional change vs time'};
exts = {'.tiff', '.eps', '.fig', '.jpg'}; ftype = {'tiff', 'eps2c', 'fig', 'jpeg'};
if (~isempty(SaveDirectory))
    s = warning('off');
    if (~isempty(SaveDirectory))
        for j = 1:4
            figure(j);
            set(j, 'InvertHardcopy', 'off');
            for k = 1:length(exts)        
                fname = fullfile(SaveDirectory, [saveName{j} exts{k}]);
                saveas(j, fname, ftype{k});
            end
        end                     
    end
    warning(s);
end

%{
figure(5); clf(5);
%dlltri_vs_etx = interp1(llet, dlltri./lltri, tri_ad_40.etxs);

[~,dlltri_vs_etx] = meanyvsx(llet, dlltri./lltri, -200:40:200);
dllexp_vs_etx = interp1(llet, dllexp./llexp, ad_exp_40.etxs);

errorbar (dlltri_vs_etx, tri_ad_40.reo_vs_toff, tri_ad_40.reo_vs_toff_eb,'bs'); hold on
errorbar (dllexp_vs_etx, ad_exp_40.reo_vs_toff, tri_ad_40.reo_vs_toff_eb, 'go');
%}


clear po;
po.onColor = [];
po.offColor = [];
po.reverse = true;

po = repmat(po, 1 ,2);
po(1).color = 'b';
po(1).marker = 's';
po(2).color = 'g';
po(2).marker = 'o';
po(1).legendEntry = 'triangle';
po(2).legendEntry = 'exponential';


existsAndDefault('savefigs', false);
whichGraphs ={'reorate_vs_time','reomag_vs_time','numhs_vs_time'};
saveDirectory = [];
sf = 5;
if (savefigs)
    saveDirectory = fullfile(sd, 'combined', 'with titles', 'with legend');
    mkdir(saveDirectory);
end
temporalNavigationFigures([tri_ad_20, ad_exp_20], po, 'showlegend', true, 'showtitle', true, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

if (savefigs)
    saveDirectory = fullfile(sd, 'combined', 'no titles', 'with legend');
    mkdir(saveDirectory);
end
temporalNavigationFigures([tri_ad_20, ad_exp_20], po, 'showlegend', true, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

if (savefigs)
    saveDirectory = fullfile(sd, 'combined', 'no titles', 'without legend');
    mkdir(saveDirectory);
end
temporalNavigationFigures([tri_ad_20, ad_exp_20], po, 'showlegend', false, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

if (savefigs)
    saveDirectory = fullfile(sd, 'triangle no background');
    mkdir(saveDirectory);
end
temporalNavigationFigures(tri_ad_20, po(1), 'showlegend', false, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

if (savefigs)
    saveDirectory = fullfile(sd, 'exponential no background');
    mkdir(saveDirectory);
end
temporalNavigationFigures(ad_exp_20, po(2), 'showlegend', false, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

po(1).offColor = 'k';
po(2).offColor = 'k';

if (savefigs)
    saveDirectory = fullfile(sd, 'triangle');
    mkdir(saveDirectory);
end
temporalNavigationFigures(tri_ad_20, po(1), 'showlegend', false, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);

if (savefigs)
    saveDirectory = fullfile(sd, 'exponential');
    mkdir(saveDirectory);
end
temporalNavigationFigures(ad_exp_20, po(2), 'showlegend', false, 'showtitle', false, 'SaveDirectory', saveDirectory, 'startFignum', sf, 'whichGraphs', whichGraphs, 'backgroundColor', 'w');
sf = sf + length(whichGraphs);
