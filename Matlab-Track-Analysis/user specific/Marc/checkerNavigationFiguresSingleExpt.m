function hset = checkerNavigationFiguresSingleExpt(ad, plot_options, varargin)
%function hset = checkerNavigationFigures(ad, plot_options, varargin)
%'HeadSwingAcceptanceHandedness_0', 'HeadSwingAcceptanceDirection_0',
%'HeadSwingAcceptanceHandedness_45','HeadSwingDirection_0,
%'HeadSwingDirection_45', 'ReorientationRateVsHeading',
%'ReorientationRateVsDistance', 'ReorientationMagVsHeading',
%'ReorientationDirVsHeading', 'PauseRateVsHeading'
% varargin:
% SaveDirectory = [];
% forprinting = false;
% showlegend = false;
% startFignum = 1;
% showtitle = false;
% whichGraphs = {};
% backgroundColor = [];
% legendLocation = 'BestOutside';

SaveDirectory = [];
forprinting = false;
showlegend = false;
startFignum = 1;
showtitle = false;
whichGraphs = {};
backgroundColor = [];
legendLocation = 'BestOutside';
varargin = assignApplicable(varargin);
fignum = startFignum -1;
plotNumber = 0;

font = 'Arial';
if (forprinting)
    fontsize = 8;
else
    fontsize = 10;
end
default_properties = get(0,'default');
font = 'Arial';
if (forprinting)
    fontsize = 8;
else
    fontsize = 10;
end

set(0,'DefaultAxesFontSize', fontsize);
set(0,'DefaultAxesFontName', font);
set(0, 'DefaultTextInterpreter', 'Latex');
set(0, 'DefaultAxesLineWidth', 2);
if (~isempty(backgroundColor))
    if (ischar(backgroundColor))
        backgroundColor = char2rgb(backgroundColor);
    end
    set(0,'DefaultFigureColor', backgroundColor);
    set(0,'DefaultTextColor', 1 - backgroundColor);
    set(0,'DefaultAxesColor',  backgroundColor);
    set(0,'DefaultAxesXColor', 1 - backgroundColor);
    set(0,'DefaultAxesYColor', 1 - backgroundColor);
   legendOptions = {'Location', legendLocation,'Interpreter','Latex', 'Color', backgroundColor, 'XColor', 1 - backgroundColor, 'YColor', 1 - backgroundColor};
else
    legendOptions = {'Location', legendLocation,'Interpreter','Latex'};
end

tolightColor = [1 0.8 0.2];
todarkColor = [0.3 0.05 0.5];
upColor = [0.1 0.7 0.3];
downColor = [0 0.1 0.9];

directions = [-180,-90, 0, 90];
colors = {todarkColor, downColor, tolightColor, upColor};



ccc = 'bgrcymk';
sss = 'sodvh>p^<';
nexp = length(ad);
for j = 1:nexp
    po(j).lineWidth = 2; %#ok<*AGROW>
    po(j).color = ccc(mod(j-1, length(ccc)) + 1);
    po(j).legendEntry = ['eset ' num2str(j)];
    po(j).marker = sss(mod(j-1, length(sss)) + 1);
    po(j).plotOptions = {};
end
if (nargin > 2 && isstruct(plot_options))
    fn = fieldnames(plot_options);
    for j = 1:length(fn)
        for k = 1:length(po)
            po(k).(fn{j}) = plot_options(k).(fn{j});
        end
    end
end

saveName =  'quadrantMap';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    r =1;
    mhs = 0.5;
     %{
    set(compass(exp(sqrt(-1) * (-1:1) * pi / 6)), 'Color', tolightColor, 'LineWidth', 4); hold on
    set(compass(exp(sqrt(-1) * (pi/2 + (-1:1) * pi / 6))), 'Color', upColor, 'LineWidth', 4);
    set(compass(exp(sqrt(-1) * (pi + (-1:1) * pi / 6))), 'Color', todarkColor, 'LineWidth', 4);
    set(compass(exp(sqrt(-1) * (3*pi/2 + (-1:1) * pi / 6))), 'Color', downColor, 'LineWidth', 4);
    %}
    t = linspace(-pi/2, pi/2,100);
    plot (r*cos(t), r*sin(t), 'k--', 'LineWidth', 3); hold on;plot (-r*cos(t), r*sin(t), 'w--', 'LineWidth', 3); 
    set(quiver([0 0 0], [0 0 0], r*cos((-1:1) * pi / 6), r*sin((-1:1) * pi / 6), 0, 'MaxHeadSize', mhs), 'Color', tolightColor, 'LineWidth', 4); 
    set(quiver([0 0 0], [0 0 0], r*cos(pi/2 + (-1:1) * pi / 6), r*sin(pi/2 + (-1:1) * pi / 6), 0, 'MaxHeadSize', mhs),'Color', upColor, 'LineWidth', 4);
    set(quiver([0 0 0], [0 0 0], r*cos(pi + (-1:1) * pi / 6), r*sin(pi + (-1:1) * pi / 6), 0, 'MaxHeadSize', mhs),'Color', todarkColor, 'LineWidth', 4);
    set(quiver([0 0 0], [0 0 0], r*cos(3*pi/2 + (-1:1) * pi / 6), r*sin(3*pi/2 + (-1:1) * pi / 6), 0, 'MaxHeadSize', mhs),'Color', downColor, 'LineWidth', 4);
    text(r*1.02, 0, '\bf{0$^{\circ}$}', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'middle', 'Color', tolightColor,'FontSize', 18,'FontWeight', 'Bold');
    text(0,r*1.02, '\bf{90$^{\circ}$}', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom', 'Color', upColor,'FontSize', 18,'FontWeight', 'Bold');
    text(-r*1.02, 0, '\bf{180$^{\circ}$}', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'middle', 'Color', todarkColor,'FontSize', 18,'FontWeight', 'Bold');
    text(0,-r*1.02, '\bf{-90$^{\circ}$}', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'Color', downColor,'FontSize', 18,'FontWeight', 'Bold');
   
    
    ax1 = gca;
    set(ax1, 'XLim', [-r*1.1, r*1.1], 'YLim', [-r*1.15 r*1.15]); axis equal;
    ax2 = cloneaxes(ax1);
    xl = get(ax2, 'Xlim');
    yl = get(ax2, 'Ylim');
    xl(2) = 0.5*(xl(1) + xl(2));
    patch ([-r -r 0 0 -r], [-r r r -r -r], [0.3 0.3 0.3]); %[xl(1) xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(2) yl(2) yl(1) yl(1)]
    xl = get(ax2, 'Xlim');
    yl = get(ax2, 'Ylim');
    xl(1) = 0.5*(xl(1) + xl(2));
    patch ([0 0 r r 0], [-r r r -r -r], [1 1 1]); %[yl(1) yl(2) yl(2) yl(1) yl(1)]
    
    set(ax2, 'Color', [0.9 0.9 0.9], 'XTick', [], 'Ytick', [])
    axes(ax1)
    set(ax1, 'Color', 'none', 'XTick', [], 'Ytick', [])
end
saveName =  'HeadSwingAcceptanceHandedness_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    if (length(ad.hs_acctoleft0) ~= 4)
        error ('you should run with head swing bin size of 90 for these plots');
    end
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx0));
    hseb =  zeros(2*nexp, length(ad.hstx0));
    
    hsbar (1:2:end, :) = ad.hs_acctoleft0; 
    hsbar (2:2:end, :) = ad.hs_acctoright0;
    hseb (1:2:end, :) = ad.hs_acctoleft0_eb; 
    hseb (2:2:end, :) = ad.hs_acctoright0_eb;
    
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    le{1} = 'to left';
    le{2} = 'to right';
    ttl = ('Head Swing Acceptances To Left/Right');
    
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl,[],[],le,[],'axis'); 
    for j = 1:4
        I = find(cosd(ad.hstx0(j)) == cosd(directions) & sind(ad.hstx0(j)) == sind(directions), 1);
        if (~isempty(I))
            set(hhh.bars(j,1), 'FaceColor', colors{I}, 'EdgeColor', colors{I});
            set(hhh.bars(j,2), 'FaceColor', 'w', 'EdgeColor',colors{I}, 'LineWidth', 3);
        end
    end
    set (hhh.bars(j+1,:), 'EdgeColor', 'k','LineWidth', 3);
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    
    if (showtitle)
        title (ttl);
        p = get(gca, 'Position');
        p(4) = 0.95 * p(4); set(gca, 'Position', p);
    end
    set(gca, 'LineWidth', 3);
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
       
end
saveName =  'HeadSwingAcceptanceDirection_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
     if (length(ad.hs_acctoleft0) ~= 4)
        error ('you should run with head swing bin size of 90 for these plots');
    end
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx0));
    hseb =  zeros(2*nexp, length(ad.hstx0));
    
    hsbar (1:2:end, :) = ad.hs_acctolight0; 
    hsbar (2:2:end, :) = ad.hs_acctodark0;
    hseb (1:2:end, :) = ad.hs_acctolight0_eb; 
    hseb (2:2:end, :) = ad.hs_acctodark0_eb;
    
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    le{1} = 'to light';
    le{2} = 'to dark';
    ttl = ('Head Swing Acceptances To Light/Dark');
    
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl,[],[],le,[],'axis'); 
    for j = 1:4
        I = find(cosd(ad.hstx0(j)) == cosd(directions) & sind(ad.hstx0(j)) == sind(directions), 1);
        if (~isempty(I))
            set(hhh.bars(j,1), 'FaceColor', colors{I}, 'EdgeColor', colors{I});
            set(hhh.bars(j,2), 'FaceColor', 'w', 'EdgeColor',colors{I}, 'LineWidth', 4);
        end
    end
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    if (showtitle)
        title (ttl);
          p = get(gca, 'Position');
        p(4) = 0.95 * p(4); set(gca, 'Position', p);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    
end

saveName =  'HeadSwingDirection_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    hsbar  = ad.firsths_dir0;
    hseb = ad.firsths_dir0_eb;
    
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
    ylbl = 'probability first head sweep is to left';
  
    ttl = ('Head Swing Initiation Bias');
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl);
    for j = 1:4
        I = find(cosd(ad.hstx0(j)) == cosd(directions) & sind(ad.hstx0(j)) == sind(directions), 1);
        if (~isempty(I))
            set(hhh.bars(j), 'FaceColor', colors{I}, 'EdgeColor', colors{I});
        end
    end
    if (showtitle)
        title (ttl);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    hold on;
    plot (get (gca, 'XLim'), [0.5 0.5], 'r--', 'LineWidth', 3);
    set(hhh.baseline, 'LineWidth', 4);
end



saveName =  'ReorientationRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    q = zeros(size(ad.txc));
    for j = 1:length(ad.txc)
        [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
        q(j) = I;
    end
    
    for j = 1:length(directions)
        inds = q == j;
        errorbar (ad.txc(inds), ad.reorate_thetabound(inds), ad.reorate_thetabound_eb(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;        
    end
    
    for j = 1:(length(ad.txc) - 1)
        xdata = interp1(ad.txc, [j j+0.5]);
        ydata = interp1(ad.reorate_thetabound, [j j+0.5]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        xdata = interp1(ad.txc, [j+0.5 j+1]);
        ydata = interp1(ad.reorate_thetabound, [j+0.5 j+1]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('reorientation rate (min$^{-1}$)');
    if(showtitle), title ('Reorientation Rate vs. Direction of Forward Movement Relative To Boundary');end 
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'ReorientationRateVsDistance';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    errorbar (10*ad.distx, ad.reorate_disttobound_todark, ad.reorate_disttobound_todark_eb, 'k-', 'Color', todarkColor, 'Marker', po.marker,  'LineWidth', po.lineWidth, po.plotOptions{:}); hold on;
    errorbar (10*ad.distx, ad.reorate_disttobound_tolight, ad.reorate_disttobound_tolight_eb, 'k-', 'Color', tolightColor, 'Marker', po.marker,  'LineWidth', po.lineWidth, po.plotOptions{:}); 
    leg{1} = 'light to dark';
    leg{2} = 'dark to light';

    xlim ([-5 5]);
    ax1 = gca;
        xlabel ('distance from boundary (mm)');
    ylabel ('reorientation rate (min$^{-1}$)');
    if (showtitle), title ('Reorientation Rate vs. Distance of Head from Boundary'); end
    if(showlegend)
        legend (leg, 'Location', 'NorthEast','Interpreter','Latex'); 
    else
        [~,I] = min(abs(ad.distx));
        quiver(10*ad.distx(I), ad.reorate_disttobound_tolight(I), 1.5, 0, 0, 'Color', tolightColor, 'LineWidth', 3, 'MaxHeadSize', 10);     
        quiver(0, 0.25*ad.reorate_disttobound_tolight(I), -1.5, 0, 0, 'Color', todarkColor, 'LineWidth', 3, 'MaxHeadSize', 10);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);

    ax2 = cloneaxes(ax1);
    xl = get(ax2, 'Xlim');
    yl = get(ax2, 'Ylim');
    xl(2) = 0.5*(xl(1) + xl(2));
    patch ([xl(1) xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(2) yl(2) yl(1) yl(1)], [0.3 0.3 0.3]);
    set(ax2, 'Color', [1 1 1], 'XTick', [], 'Ytick', [])
    axes(ax1)
    set(ax1, 'Color', 'none');
end



saveName =  'ReorientationMagVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    
    q = zeros(size(ad.txc));
    for j = 1:length(ad.txc)
        [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
        q(j) = I;
    end
    
    for j = 1:length(directions)
        inds = q == j;
        errorbar (ad.txc(inds), ad.reosize_vs_ttb(inds), ad.reosize_vs_ttb_eb(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;        
    end
    
    for j = 1:(length(ad.txc) - 1)
        xdata = interp1(ad.txc, [j j+0.5]);
        ydata = interp1(ad.reosize_vs_ttb, [j j+0.5]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        xdata = interp1(ad.txc, [j+0.5 j+1]);
        ydata = interp1(ad.reosize_vs_ttb, [j+0.5 j+1]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        
    end
    
    
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('mean reorientation magnitude (deg)');
    if (showtitle),title ('Size of Reorientation vs. Direction of Forward Movement Relative To Boundary');end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'ReorientationDirVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    q = zeros(size(ad.txc));
    for j = 1:length(ad.txc)
        [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
        q(j) = I;
    end
    
    for j = 1:length(directions)
        inds = q == j;
        errorbar (ad.txc(inds), ad.reodir_vs_ttb(inds), ad.reodir_vs_ttb_eb(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;        
    end
    
    for j = 1:(length(ad.txc) - 1)
        xdata = interp1(ad.txc, [j j+0.5]);
        ydata = interp1(ad.reodir_vs_ttb, [j j+0.5]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        xdata = interp1(ad.txc, [j+0.5 j+1]);
        ydata = interp1(ad.reodir_vs_ttb, [j+0.5 j+1]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('mean reorientation (deg)');
    if (showtitle),title ('Direction of Reorientation vs. Direction of Forward Movement Relative To Boundary');end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'PauseRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    q = zeros(size(ad.txc));
    for j = 1:length(ad.txc)
        [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
        q(j) = I;
    end
    
    for j = 1:length(directions)
        inds = q == j;
        errorbar (ad.txc(inds), ad.pauserate_thetabound(inds), ad.pauserate_thetabound(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;        
    end
    
    for j = 1:(length(ad.txc) - 1)
        xdata = interp1(ad.txc, [j j+0.5]);
        ydata = interp1(ad.pauserate_thetabound, [j j+0.5]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        xdata = interp1(ad.txc, [j+0.5 j+1]);
        ydata = interp1(ad.pauserate_thetabound, [j+0.5 j+1]);
        plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;   
        
    end
    
   
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('pausing rate (min$^{-1}$)');    
    if (showtitle), title ('Pausing Rate vs. Direction of Forward Movement Relative To Boundary');end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end


exts = {'.tiff', '.eps', '.ai', '.fig', '.jpg', '.pdf'}; ftype = {'tiff', 'eps2c', 'ai', 'fig', 'jpeg', 'pdf'};
if (~isempty(SaveDirectory))
    s = warning('off');
    for j = 1:length(hset);
        figure(hset(j).fignum);
        set(hset(j).fignum, 'InvertHardcopy', 'off');
        for k = 1:length(exts)        
            fname = fullfile(SaveDirectory, [hset(j).saveName exts{k}]);
            saveas(hset(j).fignum, fname, ftype{k});
        end
    end
    warning(s);
end
set(0,default_properties);

function nextfigure()
fignum = evalin('caller', 'fignum+1;'); assignin('caller', 'fignum', fignum);
plotNumber = evalin('caller', 'plotNumber+1;'); assignin('caller', 'plotNumber', plotNumber);

forprinting = evalin('caller', 'forprinting;');
if(forprinting)
    figureForPrinting(fignum);
else
    figure(fignum); clf(fignum);
end
