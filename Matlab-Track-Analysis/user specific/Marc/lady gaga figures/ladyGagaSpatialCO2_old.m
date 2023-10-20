%run co2analysis20101223 first
%{
if (~exist('navind_eb_small', 'var'))
    for j = 1:9
        %{
        v{j} = co2_eset(j).gatherField('vel','mean');
        sp{j} = co2_eset(j).gatherField('speed', 'mean');
        nt(j) = length(sp{j});
        np{j} = co2_eset(j).gatherField('npts');
        alpha = [np{j};np{j}]/sum(np{j});
        savg = sum(alpha(1,:).*sp{j});
        vavg = sum(alpha.*v{j}, 2);
        vvar = sum(alpha.*(v{j} - repmat(vavg, 1, length(v{j}))).^2);
        nmeas = sum(np{j}) * co2_eset(j).expt(1).track(1).dr.interpTime / (2 ( co2_eset(j).autocorr_tau));
        navindeb(:,j) = vavg / savg;
        
%        navind(:,j) = sum([np{j};np{j}].*v{j},2)/sum(np{j}.*sp{j});
 %       navindeb(:,j) = std(v{j}./[sp{j};sp{j}], 0, 2) / sqrt(nt(j));%not quite right
        %}
        v = co2_eset(j).gatherField('vel');
        s = sqrt(sum(v.^2));
        navind(:,j) = mean(v,2)/mean(s);
        nt = sum(co2_eset(1).evaluateTrackExpression('1'));
        
        navind_eb_big(:,j) = std(v,0,2)/mean(s)/sqrt(nt);
        navind_eb_small(:,j) = std(v,0,2)/mean(s)/sqrt(length(s)*co2_eset(j).expt(1).track(1).dr.interpTime/(2*co2_eset(j).autocorr_tau));

    end
end
%}
sd = 'C:\Users\Marc\Documents\figures\talk fri mar 4 2011 lady gaga\spatial navigation';
if (~exist ('sno', 'var'))
    sno.angleBinSize = 45;
    sno.minHSTheta = 20;
    sno.preferredDirection = 0;     
end
if (~exist('driftvel', 'var'))
    for j = 1:length(co2_eset)
        driftvel{j} = zeros(2, length(co2_eset(j).expt));
        for k = 1:length(co2_eset(j).expt)
            driftvel{j}(:,k) = mean(co2_eset(j).expt(k).gatherField('vel'),2);
        end
    end
end
if (false && ~exist('adlgbins', 'var'))
    snolgbins = sno;
    snolgbins.angleBinSize = 90;
    adlgbins = spatialNavigationMaggotAnalysis(ad, snolgbins);
end
ccc = 'bgrcmyw';
sss = 'sodvh>p^<';
le = {'CS CO$_2$ 0-5$\%$', 'Gr63a CO$_2$ 0-5$\%$'  ,  'CS CO$_2$ 0-1$\%$'  ,  'CS CO$_2$ 0-0.5$\%$'  ,  'CS air control'  ,  'CS no air control'  ,  'CS CO$_2$ 0-2.5$\%$'  ,  'CS Ethyl Acetate',...
    'Gr63a Ethyl Acetate'};

for j = 1:length(ad)
    po(j).lineWidth = 2;
    po(j).color = ccc(mod(j-1, length(ccc)) + 1);
    po(j).marker = sss(mod(j-1, length(sss)) + 1);
    po(j).plotOptions = {};
    po(j).legendEntry = le{j};
end
po(8).color = 'r';
existsAndDefault('savefiles', false);
font = 'Arial';
fontsize = 10;

set(0,'DefaultAxesFontSize', fontsize);
set(0,'DefaultAxesFontName', font);
set(0,'DefaultTextInterpreter', 'Latex');

esetinds = [1 8 2 5];
thinfactor = [1 1 1 1];% [0.5 0.3 0.8 1];
if (exist ('co2_eset', 'var'))
    figure(51);
    clf();
    for j = 1:4
        subplot(2,2,j);
        cla();
        hold on;
        co2_eset(esetinds(j)).evaluateTrackExpression(['if (rand(1) < ' num2str(thinfactor(j)) ' && (track.npts > 2000 || any(sum(track.getDerivedQuantity(''displacement'').^2) > 4.5^2)) && min(track.getDerivedQuantity(''eti'')) < 600), track.plotPath(''displacement'',''k-'',''Color'',[0.3 0.3 0.3]); end;']);
        plot (0,0,'r.','MarkerSize', 20);
        %[~,ml] = co2_eset(esetinds(j)).meanField2vsField1('eti', 'sloc', 120:10:1020);
        %plot (ml(1,:)-ml(1,1), ml(2,:) - ml(2,1), 'k-', 'LineWidth', 4);
        title (po(esetinds(j)).legendEntry);
        axis equal;
        axis([-12 12 -6 6]);
        set(gca, 'XTick', [],'YTick', [], 'box', 'on','LineWidth',2)
        emsmallen(gca, 'FontSize', fontsize, 'Font', font);     

    end
    if (false && savefiles)
        saveas(gcf, fullfile(sd, 'tracks.fig'), 'fig');
        saveas(gcf, fullfile(sd, 'tracks.eps'), 'eps2c');
        saveas(gcf, fullfile(sd, 'tracks.jpg'), 'fig');
        saveas(gcf, fullfile(sd, 'tracks.tif'), 'tiff');
    end
end
%{
figure(52);
clf();
order = [1 7 3 4 5 6 2 8 9];
% barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
 
barhandles = barweb(navind(:,order)', navind_eb_small(:,order)', [], {po(order).legendEntry}, [], [], 'Navigation Index');
yl = get(barhandles.ax, 'YLim');
set(barhandles.ax, 'YLim', [-max(abs(yl)) max(abs(yl))], 'box', 'on','FontSize', fontsize);
set(barhandles.ax, 'Position', get(barhandles.ax,'Position') + [0 0 0.05 0]);
legend ('Parallel to Gradient', 'Perpendicular to Gradient', 'Location', 'Best');
title ('Navigation of Spatial Carbon Dioxide Gradients');
emsmallen(barhandles.ax, 'FontSize', fontsize, 'Font', font);
th = rotateticklabel(barhandles.ax, -20,false);
set(th, 'FontSize', fontsize-2);

%}
backgroundColor = [0 0 0];

set(0,'DefaultFigureColor', backgroundColor);
set(0,'DefaultTextColor', 1 - backgroundColor);
set(0,'DefaultAxesColor',  backgroundColor);
set(0,'DefaultAxesXColor', 1 - backgroundColor);
set(0,'DefaultAxesYColor', 1 - backgroundColor);

figure(53);
clf();
esetinds = [1 7 8 2 5 6];
for j = 1:6
    subplot(2,3,j);
    cla();
  %  set (gcf, 'Color', backgroundColor);
   % set (gca, 'Color', backgroundColor, 'XColor', 1-backgroundColor, 'YColor', 1 - backgroundColor);
    h = polar (1/8,0.6, 'k'); hold on;
    %set(get(h,'Parent'), 'Color', backgroundColor, 'XColor', 1-backgroundColor, 'YColor', 1 - backgroundColor);
    
    if (po(esetinds(j)).color == 'k')
        c = 'w';
    else
        c =  po(esetinds(j)).color;
    end
    set(compass(60*driftvel{esetinds(j)}(1,:), 60*driftvel{esetinds(j)}(2,:),c),'LineWidth', 2);
    hold off
    title ([po(esetinds(j)).legendEntry ' - ' num2str(length(driftvel{esetinds(j)})) ' expts']);
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    
end

if (savefiles)
    set(gcf, 'InvertHardcopy', 'off');
    saveas(gcf, fullfile(sd, 'drift vel.fig'), 'fig');
    saveas(gcf, fullfile(sd, 'drift vel.eps'), 'eps2c');
    saveas(gcf, fullfile(sd, 'drift vel.tif'), 'tiff');
end
%esetinds = [1 2 5];
%spatialNavigationMaggotFigures(ad(esetinds), sno, po(esetinds));

esetinds = [1 7 3 4 5 6 2];
whichGraphs = {'NavigationIndex', 'StrategicIndices'};
if (savefiles)
    savedir = fullfile (sd, 'all co2');
    mkdir (savedir);
else
    savedir = [];
end
spatialNavigationMaggotFigures(ad(esetinds), sno, po(esetinds), 'whichGraphs', whichGraphs, 'SaveDirectory', savedir, 'showlegend', true, 'backgroundColor', backgroundColor);

fignum = length(whichGraphs) + 1;
esetinds = [1 2 5 6];
whichGraphs = {'NavigationIndex', 'DirectionHistogram',  'RunStartHistogram', 'InstantaneousDeltaThetaVsTheta', 'SpeedVsDirection', 'ReoDirVsHeading','ReoMagVsHeading','HeadSwingAcceptanceDirection'};
if (savefiles)
    savedir = fullfile (sd, 'co2 controls');
    mkdir (savedir);
else
    savedir = [];
end
spatialNavigationMaggotFigures(ad(esetinds), sno, po(esetinds), 'whichGraphs', whichGraphs, 'SaveDirectory', savedir, 'startFigNum', fignum, 'showlegend', true, 'backgroundColor', backgroundColor);


fignum = fignum + length(whichGraphs);
esetinds = [1 8];
whichGraphs = {'NavigationIndex', 'ReorientationRateVsHeading','StrategicIndices', 'DirectionHistogram',  'RunStartHistogram', 'InstantaneousDeltaThetaVsTheta', 'SpeedVsDirection', 'ReoDirVsHeading','ReoMagVsHeading','HeadSwingAcceptanceDirection'};
if (savefiles)
    savedir = fullfile (sd, 'ethyl acetate vs co2');
    mkdir (savedir);
else
    savedir = [];
end
spatialNavigationMaggotFigures(ad(esetinds), sno, po(esetinds), 'whichGraphs', whichGraphs, 'SaveDirectory', savedir, 'startFigNum', fignum, 'showlegend', true, 'backgroundColor', backgroundColor);

set(0,'DefaultFigureColor', 'remove');
set(0,'DefaultTextColor', 'remove');
set(0,'DefaultAxesColor',  'remove');
set(0,'DefaultAxesXColor', 'remove');
set(0,'DefaultAxesYColor', 'remove');

