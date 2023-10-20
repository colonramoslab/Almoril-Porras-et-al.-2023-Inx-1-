function figurehandles = spatialNavigationMaggotFigures (esets, spatial_navigation_options, plot_options, varargin)
%function spatialNavigationMaggotFigures (esets, spatial_navigation_options, plot_options, varargin)
%
% spatial_navigation_options -- 1 set of options for all esets
% plot_options -- each eset gets its own options
% plot_options.
%    lineWidth -- width of lines
%    color -- color for line and marker
%    legendEntry -- what to put in legend
%    marker -- marker 
%    plotOptions -- additional options to pass to plotting function
%
%   varargin: 'SaveDirectory', directory to save images in (if empty,
%   nothing saved
%   forprinting = false;
%   showlegend = false;
%   vectorgraphics = whether to enforce vector based graphics (painters);  defaults to
%   true unless we are not saving figures or showing legends, in which case
%   we use openGL renderer, allowing transparency
%   whichGraphs = {}; all if empty - choose from {'DirectionHistogram',
%   'RunStartHistogram','ReorientationRateVsHeading',
%   'InstantaneousDeltaThetaVsTheta','SpeedVsDirection',
%   'ReoDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness',
%   'HeadSwingAcceptanceDirection', 'HeadSwingDirection',
%   'FirstHeadSwingHandedness', 'FirstHeadSwingMeanDir', 'NavigationIndex',
%   'StrategicIndices','FirstHeadSwingBias'
%   startFigNum = 1; 

SaveDirectory = [];
forprinting = false;
showlegend = false;
vectorgraphics = [];
whichGraphs = {'DirectionHistogram', 'ReorientationRateVsHeading','SpeedVsDirection','ReoDirVsHeading','ReoMagVsHeading', 'HeadSwingAcceptanceHandedness','NavigationIndex', 'StrategicIndices'};
%ignoring
%   'RunStartHistogram',
%   'InstantaneousDeltaThetaVsTheta',    
%   'HeadSwingAcceptanceDirection', 'HeadSwingDirection',
%   'FirstHeadSwingHandedness', 'FirstHeadSwingMeanDir', ,'FirstHeadSwingBias'
startFigNum = 1;
backgroundColor = [];
legendLocation = 'BestOutside';
varargin = assignApplicable(varargin);
if (isempty(vectorgraphics))
    vectorgraphics = ~(isempty(SaveDirectory) && ~showlegend);
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
sno.angleBinSize = 30; % in degrees
sno.preferredDirection = 0;
sno.hsBinSize = 90;

set(0, 'DefaultTextInterpreter', 'Latex');
%set(0, 'DefaultLegendInterpreter', 'Latex');
if (~isempty(backgroundColor))
    if (ischar(backgroundColor))
        backgroundColor = char2rgb(backgroundColor);
    end
    set(0,'DefaultFigureColor', backgroundColor);
    set(0,'DefaultTextColor', 1 - backgroundColor);
    set(0,'DefaultAxesColor',  backgroundColor);
    set(0,'DefaultAxesXColor', 1 - backgroundColor);
    set(0,'DefaultAxesYColor', 1 - backgroundColor);
   % set(0,'DefaultLegendColor',  1 - backgroundColor);
   % set(0,'DefaultLegendXColor', 1 - backgroundColor);
   % set(0,'DefaultLegendYColor', 1 - backgroundColor);
   legendOptions = {'Location', legendLocation,'Interpreter','Latex', 'Color', backgroundColor, 'XColor', 1 - backgroundColor, 'YColor', 1 - backgroundColor};
else
    legendOptions = {'Location', legendLocation,'Interpreter','Latex'};
end

if (nargin == 0)
    figurehandles = sno;
    return;
end

if (nargin > 1 && isstruct(spatial_navigation_options))
    fn = fieldnames(spatial_navigation_options);
    for j = 1:length(fn)
        sno.(fn{j}) = spatial_navigation_options.(fn{j});
    end
end

if (isa(esets, 'ExperimentSet'))
    ad = spatialNavigationMaggotAnalysis(esets, sno);
else
    ad = esets;
end


%{
fields = fieldnames(analyzedData);

for j = 1:length(fields)
    eval([fields{j} ' = analyzedData.' fields{j}  ';']);
end
nexp = size(reohist, 1);
%}
ccc = 'bgrcymk';
sss = 'sodvh>p^<';
nexp = length(ad);
for j = 1:nexp
    po(j).lineWidth = 2;
    po(j).color = ccc(mod(j-1, length(ccc)) + 1);
    po(j).legendEntry = ['eset ' num2str(j)];
    po(j).marker = sss(mod(j-1, length(sss)) + 1);
    po(j).plotOptions = {};
    po(j).shadedErrorRegion = false;
    po(j).useGauss = false;
end

if (nargin > 2 && isstruct(plot_options))
    fn = fieldnames(plot_options);
    for j = 1:length(fn)
        for k = 1:length(po)
            po(k).(fn{j}) = plot_options(k).(fn{j});
        end
    end
end
for j = 1:length(po)
    if (ischar(po(j).color))
        po(j).color = char2rgb(po(j).color);
    end
end
% 
 if (showlegend && ~vectorgraphics)
     disp ('warning:  cannot show legends with opengl.  using painters renderer, which does not have transparency');
 end
%     [po.shadedErrorRegion] = deal(false);
% end
%{
tx = ad(1).tx;
txc = [tx tx(end)+sno.angleBinSize]; %#ok<*COLND>
txe = txc - sno.angleBinSize/2;
%}
fignum = startFigNum-1;
plotNumber = 0;
set(0, 'DefaultTextInterpreter', 'Latex');
nexp = length(ad);
%fignum = fignum + 1; figure(fignum); clf(fignum);
saveName = 'DirectionHistogram';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
  
    
    
    errorplot(ad,'txc', 'thetahist', 'thetahist_eb', po, 1./sum(reshape([ad.thetahist],[],length(ad))), true, showlegend);
    
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('relative frequency');
    title ('Instantaneous Direction of Forward Movement');
    
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    
end

saveName =  'ReorientationRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    errorplot(ad,'txc', 'reohist', 'reohist_eb', po,1, true, showlegend);
    %{
    for j = 1:nexp
        errorbar (ad(j).txc, ad(j).reohist([1:end 1]),  ad(j).reohist_eb([1:end 1]), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
    end
    %}
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('reorientation rate (min$^{-1}$)');
    title ('Reorientation Rate vs. Direction of Forward Movement');
    %if(showlegend), legend ({po.legendEntry}, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'RunStartHistogram';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    errorplot(ad,'txc', 'runStartDirectionHist', 'runStartDirectionHist_eb', po,1, true, showlegend);
    
% 
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).runStartDirectionHist([1:end 1]),ad(j).runStartDirectionHist_eb([1:end 1]),'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('start direction of run (degrees)');
    ylabel ('relative frequency');
    title ('Heading at beginning of runs');
%     if(showlegend), legend ({po.legendEntry},legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end
saveName =  'InstantaneousDeltaThetaVsTheta';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    yfield = 'instantaneousdthetavstheta';
    ymult = 60*180/pi;
    
    if (po(1).useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'txc';
        iscirc = true;
    end
    
    errorplot(ad,xfield, yfield, [yfield '_eb'], po,ymult, iscirc, showlegend);
    
%     for j = 1:nexp
%         errorbar (ad(j).txc, 60*rad2deg(ad(j).instantaneousdthetavstheta([1:end 1])), 60*rad2deg(ad(j).instantaneousdthetavstheta_eb([1:end 1])),'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    title ('Heading Change vs. Direction (in runs)');
    xlabel ('heading (degrees)');
    ylabel ('mean heading change (degrees/min)');
%     if (showlegend), legend ({po.legendEntry},legendOptions{:});end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'SpeedVsDirection';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    yfield = 'speedVsDir';
    ymult = 1;
    
    if (po(1).useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'txc';
        iscirc = true;
    end
    
    errorplot(ad,xfield, yfield, [yfield '_eb'], po,ymult, iscirc, showlegend);
   % errorplot(ad,'txc', 'speedVsDir', 'speedVsDir_eb', po,1, true, showlegend);
    

%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).speedVsDir([1:end 1]), ad(j).speedVsDir_eb([1:end 1]), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('heading (degrees)');
    ylabel ('speed');
    title ('Speed vs. Direction (in runs)');
%     if(showlegend), legend ({po.legendEntry}, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end
saveName =  'ReoDirVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    yfield = 'reoDir';
    ymult = 180/pi;
    
    if (po(1).useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'reotxc';
        iscirc = true;
    end
    
    errorplot(ad,xfield, yfield, [yfield '_eb'], po,ymult, iscirc, showlegend);
%    errorplot(ad,'reotxc', 'reoDir', 'reoDir_eb', po,180/pi, true, showlegend);
    
%     
%     for j = 1:nexp
%         errorbar (ad(j).reotxc, rad2deg(ad(j).reoDir([1:end 1])), rad2deg(ad(j).reoDir_eb([1:end 1])), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('previous heading (degrees)');
    ylabel ('mean reorientation direction');
    title ('Mean reorientations vs. previous heading');
    if(showlegend), legend ({po.legendEntry}, legendOptions{:}); end

%     emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end
saveName =  'ReoMagVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    
    yfield = 'reoMag';
    ymult = (180/pi)^2;
    
    if (po(1).useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'reotxc';
        iscirc = true;
    end
    
    errorplot(ad,xfield, yfield, [yfield '_eb'], po,ymult, iscirc, showlegend);
%        errorplot(ad,'reotxc', 'reoMag', 'reoMag_eb', po,(180/pi).^2, true, showlegend);

% 
%     for j = 1:nexp
%         errorbar (ad(j).reotxc, rad2deg(ad(j).reoMag([1:end 1])), rad2deg(ad(j).reoStd([1:end 1])), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('previous heading (degrees)');
    hy = ylabel ('$<\Delta\theta^2>$');
    title ('Mean-Square of reorientations vs. previous heading');
   
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    yt = get(gca, 'YTick');
    yt = unique((5*round((sqrt(yt)/5))));
    pos = get(hy, 'Position');
    set (gca, 'YTick', []);
    xl = get(gca, 'XLim');
    
    for j = 1:length(yt)
        ht(j) = text(xl(1), yt(j).^2, ['$(' num2str(yt(j)) '^\circ)^2$']);
    end
    set(ht, 'FontSize', fontsize, 'FontName', font, 'HorizontalAlignment', 'Right');
    set(hy, 'Position', pos);
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'HeadSwingAcceptanceHandedness';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;


    hsbar = zeros(2*nexp, length(ad(1).hstx));
    for j = 1:nexp
        hsbar (2*j-1, :) = ad(j).headSwingAcceptanceRateLeft;
        hsbar (2*j, :) = ad(j).headSwingAcceptanceRateRight;
    end
    hseb =  zeros(2*nexp, length(ad(1).hstx));
    for j = 1:nexp
        hseb (2*j-1, :) = ad(j).headSwingAcceptanceRateLeft_eb;
        hseb (2*j, :) = ad(j).headSwingAcceptanceRateRight_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' to left'];
        le{(2*j)} = [po(j).legendEntry ' to right'];
    end
    ttl = ('Head Swing Acceptances By Bend Direction');
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    for j = 1:nexp
        set(hhh.bars(2*j - 1), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
        set(hhh.bars(2*j), 'FaceColor', 'w', 'EdgeColor', po(j).color, 'LineWidth', 4);
    end

    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
     pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

saveName =  'HeadSwingAcceptanceDirection';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;


    hsbar = zeros(nexp, 2);
    hsbar(:, 1) = [ad.headSwingAcceptanceRateTowards];
    hsbar(:, 2) = [ad.headSwingAcceptanceRateAway];
    hseb =  hsbar;
    hseb(:, 1) = [ad.headSwingAcceptanceRateTowards_eb];
    hseb(:, 2) = [ad.headSwingAcceptanceRateAway_eb];
    gnames = {po.legendEntry};
    ylbl = 'probability of accepting head swing';
    ttl = ('Head Swing Acceptances By Bend Direction');
    bardat = hsbar;
    ebdat = hseb;
    
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
    th = rotateticklabel(hhh.ax, 20);
    set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['toward ' num2str(mod(sno.preferredDirection + 90, 360)) '$^{\circ}$']), legendOptions{:});
 %   set(hhh.ax, 'YLim', [-1 1]);
 
     if (isempty(backgroundColor))
         bc = [1 1 1];
     else
         bc = backgroundColor;
     end
     
    for j = 1:length(ad)
        if ischar(po(j).color)
            cc = char2rgb(po(j).color);
        else
            cc = po(j).color;
        end
 %       po(j).color
%        0.2*po(j).color + 0.8*bc
        set(hhh.bars(j,1), 'FaceColor', 0.7*cc + 0.3*bc, 'EdgeColor', po(j).color);
        set(hhh.bars(j,2), 'FaceColor', 0.2*cc + 0.8*bc, 'EdgeColor', po(j).color,'LineWidth',1.5);
    end
    
    %hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],[],[]);
    %th = rotateticklabel(hhh.ax, 20);
    %set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['away from ' num2str(sno.preferredDirection) '$^{\circ}$']),legendOptions{:});
    set(hhh.ax, 'YLim', [0 1]);
     pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    %{
    if (isempty(backgroundColor))
        set (hhh.bars(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', [0.7 0.7 0.7]);
        set (hhh.bars(2), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 3);
    else
        set (hhh.bars(1), 'FaceColor', 0.7 - 0.7*backgroundColor, 'EdgeColor', 0.3 - 0.3*backgroundColor);
        set (hhh.bars(2), 'FaceColor', 0.3 - 0.3*backgroundColor, 'EdgeColor', 0.7 - 0.7*backgroundColor, 'LineWidth', 3);
    end
    %}
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

saveName =  'HeadSwingDirection';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;


    hsbar = zeros(2*nexp, length(ad(j).hstx));
    for j = 1:nexp
        hsbar (2*j-1, :) = rad2deg(ad(j).meanAcceptedHeadSwingDir);
        hsbar (2*j, :) = rad2deg(ad(j).meanRejectedHeadSwingDir);
    end
    hseb =  zeros(2*nexp, length(ad(j).hstx));
    for j = 1:nexp
        hseb (2*j-1, :) = rad2deg(ad(j).meanAcceptedHeadSwingDir_eb);
        hseb (2*j, :) = rad2deg(ad(j).meanRejectedHeadSwingDir_eb);
    end
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    ylbl = 'mean headswing angle';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' accepted'];
        le{(2*j)} = [po(j).legendEntry ' rejected'];
    end
    ttl = ('Mean Headswing Angle vs. Tail Direction');
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');

    for j = 1:nexp
        set(hhh.bars(2*j - 1), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
        set(hhh.bars(2*j), 'FaceColor', 'w', 'EdgeColor', po(j).color, 'LineWidth', 4);
    end
    
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
     pos = get(hhh.ax, 'OuterPosition');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'OuterPosition', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

saveName =  'FirstHeadSwingHandedness';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    ttl = 'First Head Swing Bias vs. Previous Direction';
    ylbl = 'probability first head sweep is to left';
    le = {po(1:nexp).legendEntry};
    %{
    bardat = firstHSDir';
    ebdat = firstHSDir_eb';
    %}
    bardat = reshape([ad.firstHSDir], [], length(ad));
    ebdat = reshape([ad.firstHSDir_eb], [], length(ad));

    %Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
    if (nexp > 1)
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    else
        gnames
        ttl
        ylbl
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
      %  set(hhh.ax, 'XTick', 1:length(bardat), 'XTickLabel', gnames);
    end
    for j = 1:nexp
        set(hhh.bars(j), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
    end
    %xlabel ('previous direction');
    %ylabel ('probability first head sweep is to left');
    %title ('First Head Sweep Directional Bias');
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (nexp > 1)
        set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    end
    hold on
    xl = get(gca, 'Xlim');
    plot (xl, [0.5 0.5], 'r--','LineWidth', 3);
    hold off
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

saveName =  'FirstHeadSwingBias';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    ttl = 'First Head Swing Bias vs. Previous Direction';
    ylbl = 'probability first hs in direction';
   
    %{
    bardat = firstHSDir';
    ebdat = firstHSDir_eb';
    %}
    numtowards = zeros(length(ad), 1);
    numaway = numtowards;
    for j = 1:length(ad)
    %firsths_all_dist{3} = hist(hsmt(leftofpref), deg2rad(hsdtx));
    %firsths_all_dist{4} = hist(hsmt(rightofpref), deg2rad(hsdtx));

        numtowards(j) = sum(ad(j).firsths_all_dist{3}(ad(j).hsdtx < 0)) + sum(ad(j).firsths_all_dist{4}(ad(j).hsdtx > 0));
        numaway(j) = sum(ad(j).firsths_all_dist{4}(ad(j).hsdtx < 0)) + sum(ad(j).firsths_all_dist{3}(ad(j).hsdtx > 0));
    end
    numtotal = numtowards + numaway;
    hsbar = [numtowards./numtotal numaway./numtotal];
    hseb = sqrt(hsbar.*(1-hsbar))./sqrt([numtotal numtotal]);
    gnames = {po.legendEntry};
    
    bardat = hsbar;
    ebdat = hseb;
    
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
    th = rotateticklabel(hhh.ax, 20);
    set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['toward ' num2str(mod(sno.preferredDirection + 90, 360)) '$^{\circ}$']), legendOptions{:});
 %   set(hhh.ax, 'YLim', [-1 1]);
 
     if (isempty(backgroundColor))
         bc = [1 1 1];
     else
         bc = backgroundColor;
     end
     
    for j = 1:length(ad)
        if ischar(po(j).color)
            cc = char2rgb(po(j).color);
        else
            cc = po(j).color;
        end
 %       po(j).color
%        0.2*po(j).color + 0.8*bc
        set(hhh.bars(j,1), 'FaceColor', 0.7*cc + 0.3*bc, 'EdgeColor', po(j).color);
        set(hhh.bars(j,2), 'FaceColor', 0.2*cc + 0.8*bc, 'EdgeColor', po(j).color,'LineWidth',1.5);
    end
    
    %hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],[],[]);
    %th = rotateticklabel(hhh.ax, 20);
    %set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    hold (hhh.ax, 'on');
    plot (hhh.ax, get(hhh.ax, 'XLim'), [0.5 0.5], 'r--', 'LineWidth', 2);
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['away from ' num2str(sno.preferredDirection) '$^{\circ}$']),legendOptions{:});
    set(hhh.ax, 'YLim', [0.3 0.7]);
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
    
end


saveName =  'FirstHeadSwingMeanDir';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    ttl = 'First Head Swing Bias vs. Previous Direction';
    ylbl = 'average direction of first head sweep';
    le = {po(1:nexp).legendEntry};
    %{
    bardat = firstHSDir';
    ebdat = firstHSDir_eb';
    %}
    bardat = double(rad2deg(reshape([ad.firstHSMeanDir], [], length(ad))));
    ebdat = double (rad2deg(reshape([ad.firstHSMeanDir_eb], [], length(ad))));

    %Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
    if (nexp > 1)
        %size(bardat)
        %size(ebdat)
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    else
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
      %  set(hhh.ax, 'XTick', 1:length(bardat), 'XTickLabel', gnames);
    end
    for j = 1:nexp
        set(hhh.bars(j), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
    end
    %legend ({po(1:nexp).legendEntry}, 'Location',
    %'NorthEast','Interpreter','Latex');
    %xlabel ('previous direction');
    %ylabel ('probability first head sweep is to left');
    %title ('First Head Sweep Directional Bias');
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (nexp > 1)
        set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    end
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end


saveName =  'FirstHeadSwingMeanDir';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    ttl = 'First Head Swing Bias vs. Previous Direction';
    ylbl = 'average direction of first head sweep';
    le = {po(1:nexp).legendEntry};
    %{
    bardat = firstHSDir';
    ebdat = firstHSDir_eb';
    %}
    bardat = double(rad2deg(reshape([ad.firstHSMeanDir], [], length(ad))));
    ebdat = double (rad2deg(reshape([ad.firstHSMeanDir_eb], [], length(ad))));

    %Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
    if (nexp > 1)
        %size(bardat)
        %size(ebdat)
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    else
        hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
      %  set(hhh.ax, 'XTick', 1:length(bardat), 'XTickLabel', gnames);
    end
    for j = 1:nexp
        set(hhh.bars(j), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
    end
    %legend ({po(1:nexp).legendEntry}, 'Location', 'NorthEast','Interpreter','Latex');
    
    %xlabel ('previous direction');
    %ylabel ('probability first head sweep is to left');
    %title ('First Head Sweep Directional Bias');
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (nexp > 1)
        set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    end
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end


saveName =  'NavigationIndex';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;


    hsbar = [ad.navind]';
    hseb = [ad.navind_eb]';
    gnames = {po.legendEntry};
    ylbl = '$\bar{v}/\bar{s}$';
    ttl = ('Navigation Indices');
    bardat = hsbar;
    ebdat = hseb;
%    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl);
    th = rotateticklabel(hhh.ax, 20);
    set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['toward ' num2str(mod(sno.preferredDirection + 90, 360)) '$^{\circ}$']), legendOptions{:});
 %   set(hhh.ax, 'YLim', [-1 1]);
 
     if (isempty(backgroundColor))
         bc = [1 1 1];
     else
         bc = backgroundColor;
     end
    for j = 1:length(ad)
        if ischar(po(j).color)
            cc = char2rgb(po(j).color);
        else
            cc = po(j).color;
        end
 %       po(j).color
%        0.2*po(j).color + 0.8*bc
        set(hhh.bars(j,1), 'FaceColor', 0.7*cc + 0.3*bc, 'EdgeColor', po(j).color);
        set(hhh.bars(j,2), 'FaceColor', 0.2*cc + 0.8*bc, 'EdgeColor', po(j).color,'LineWidth',1.5);
    end
 %{
    if (isempty(backgroundColor))
        set (hhh.bars(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', [0.7 0.7 0.7]);
        set (hhh.bars(2), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 3);
    else
        set (hhh.bars(1), 'FaceColor', 0.7 - 0.7*backgroundColor, 'EdgeColor', 0.3 - 0.3*backgroundColor);
        set (hhh.bars(2), 'FaceColor', 0.3 - 0.3*backgroundColor, 'EdgeColor', 0.7 - 0.7*backgroundColor, 'LineWidth', 3);
    end
    %}
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
    p = get(hhh.ax, 'Position');
    r = 0.15;
    p(1) = p(1) + r * p(3);
    p(3) = p(3) * (1-r);
    p(2) = p(2) + r * p(4);
    p(4) = (1-r) * p(4);
    set(hhh.ax, 'Position', p);
end



saveName =  'StrategicIndices';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    

    hsbar = [[ad.reoRateIndex];[ad.reoMagIndex];[ad.reoDirIndex];[ad.headSwingAcceptanceRateIndex]];
    hseb = [[ad.reoRateIndex_eb];[ad.reoMagIndex_eb];[ad.reoDirIndex_eb];[ad.headSwingAcceptanceRateIndex_eb]];
    gnames = {{'reorientation','rate'}, {'reorientation','magnitude'},{'reorientation','direction'} {'head swing', 'acceptance'}};
    le = {po.legendEntry};
    ylbl = 'index';
    ttl = ('Navigation Indices');
    bardat = hsbar;
    ebdat = hseb;
   % hhh = barweb(bardat, ebdat, 0.8, [],ttl,  [],ylbl,[],[],le,[],'axis');
    if (nexp > 1)
        hhh = barweb(bardat, ebdat, 0.8, [],ttl,  [],ylbl,[],[],le,[],'axis');
    else
        hhh = barweb(bardat', ebdat', 0.8, [],ttl,  [],ylbl);
        set(hhh.bars,  'FaceColor', po.color, 'EdgeColor', po.color);
      %  set(hhh.ax, 'XTick', 1:length(bardat), 'XTickLabel', gnames);
    end
    if (nexp > 1)
        set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    end
  
    xt = get(hhh.ax, 'XTick');
    yl = repmat(max(get(hhh.ax, 'YLim')), size(xt));
    for j = 1:length(xt)
        text (xt(j), yl(j), gnames{j}, 'FontSize', fontsize, 'FontName', font, 'VerticalAlignment', 'top','HorizontalAlignment','center');
    end
    for j = 1:nexp
        set(hhh.bars(j), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
    end
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    p = get(hhh.ax, 'Position');
    r = 0.15;
    p(1) = p(1) + r * p(3);
    p(3) = p(3) * (1-r);
    p(2) = p(2) + r * p(4);
    p(4) = (1-r) * p(4);
    set(hhh.ax, 'Position', p);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

% 
% if (~isempty(SaveDirectory))
%     for j = 1:length(hset);
%         figure(hset(j).fignum);
%         set(hset(j).fignum, 'InvertHardcopy', 'off');
%         fname = fullfile(SaveDirectory, [hset(j).saveName '.tiff']);
%         saveas(hset(j).fignum, fname, 'tiff');
%         fname = fullfile(SaveDirectory, [hset(j).saveName '.eps']);
%         saveas(hset(j).fignum, fname, 'eps2c');
%         fname = fullfile(SaveDirectory, [hset(j).saveName '.fig']);
%         saveas(hset(j).fignum, fname, 'fig');
%         
%     end
% end
if (vectorgraphics)
    exts = {'.fig', '.pdf', '.ai', '.eps'}; 
    ftype = {'fig', 'pdf', 'ill', 'eps2c'}; 
else
    exts = {'.fig', '.tiff', '.jpg'}; 
    ftype = {'fig','tiff', 'jpeg'}; 
end
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

figurehandles = hset;
set(0, default_properties);

function hh = errorplot (ad, xfield, yfield, ebfield, po, ymult, iscirc, showlegend)
    existsAndDefault('ymult', 1);
    existsAndDefault('iscirc', false);
    existsAndDefault('showlegend', false);
    nexp = length(ad);
    vg = showlegend || evalin('caller', 'vectorgraphics');
    sd = evalin('caller', '~isempty(SaveDirectory)');
    
    if (vg)
        set(gcf, 'renderer', 'painters');
    else
        set(gcf, 'renderer', 'openGL');
    end
    
    if (po(1).shadedErrorRegion)
            
        xd = {ad.(xfield)};
        for j = 1:nexp
            if (iscirc)
                inds = [1:length(ad(j).(yfield)) 1];
            else
                inds = 1:length(ad(j).(yfield));
            end
            if (length(ymult) >= j)
                ym = ymult(j);
            end
            yd{j} = ad(j).(yfield)(inds)*ym;
            eb{j} =  ad(j).(ebfield)(inds)*ym;
        end
        
        
       
        
        hh = shadedErrorPlot (xd,yd,eb,[]); 
        for j = 1:nexp
            if (~vg)
                set(hh(j), 'FaceColor', po(j).color, 'FaceAlpha', 0.25);
            else
                if (sd)
                    set(hh(j), 'FaceColor', po(j).color);
                end
            end
            set(hh(j + nexp), 'Color', po(j).color, 'Marker', po(j).marker, 'LineWidth', po(j).lineWidth, po(j).plotOptions{:});
        end
        if (showlegend)
            legendOptions = evalin('caller', 'legendOptions');
            legend (hh((nexp+1):end),{po.legendEntry}, legendOptions{:}); 
        end
        
        
        
    else
        for j = 1:nexp
            if (iscirc)
                inds = [1:length(ad(j).(yfield)) 1];
            else
                inds = 1:length(ad(j).(yfield));
            end
            if (length(ymult) >= j)
                ym = ymult(j);
            end
            hh = errorbar (ad(j).(xfield),  ad(j).(yfield)(inds)*ym, ad(j).(ebfield)(inds)*ym, 'Color', po(j).color, 'Marker', po(j).marker, 'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
        end
        if (showlegend)
            legendOptions = evalin('caller', 'legendOptions');
            legend ({po.legendEntry}, legendOptions{:}); 
        end
    end

    
function nextfigure()
    fignum = evalin('caller', 'fignum+1;'); assignin('caller', 'fignum', fignum);
    plotNumber = evalin('caller', 'plotNumber+1;'); assignin('caller', 'plotNumber', plotNumber);
    bc = evalin('caller', 'backgroundColor');
   
    forprinting = evalin('caller', 'forprinting;');
    if(forprinting)
        figureForPrinting(fignum);
    else
        figure(fignum); clf(fignum);
    end
    if (~isempty(bc))
        set(fignum, 'Color', bc);
    end
    %set(fignum, 'Color', 'w');