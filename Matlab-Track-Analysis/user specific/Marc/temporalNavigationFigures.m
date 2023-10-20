function hset = temporalNavigationFigures (ad, plot_options, varargin)
%function hset = temporalNavigationFigures (ad, plot_options, varargin)
%
%optional args:
% SaveDirectory = [];
% forprinting = false;
% showlegend = false;
% startFignum = 1;
% whichGraphs = {}; 'reorate_vs_time','speed_vs_time','reomag_vs_time','numhs_vs_time','val_vs_time'
% legendLocation = 'BestOutside';
% showtitle = true;

SaveDirectory = [];
forprinting = false;
showlegend = false;
startFignum = 1;
whichGraphs = {};
legendLocation = 'BestOutside';
showtitle = true;
varargin = assignApplicable(varargin);

plotNumber = 0;
fignum = startFignum -1;

font = 'Arial';
if (forprinting)
    fontsize = 8;
else
    fontsize = 10;
end

set(0,'DefaultAxesFontSize', fontsize);
set(0,'DefaultAxesFontName', font);
set(0, 'DefaultTextInterpreter', 'Latex');

ccc = 'bgrcymk';
sss = 'sodvh>p^<';
if (nargin == 0)
    nexp = 1;
else
    nexp = length(ad);
end

for j = 1:nexp
    po(j).lineWidth = 2;
    po(j).color = ccc(mod(j-1, length(ccc)) + 1);
    po(j).legendEntry = ['eset ' num2str(j)];
    po(j).marker = sss(mod(j-1, length(sss)) + 1);
    po(j).plotOptions = {};
    po(j).backgroundColor = [];
    po(j).onColor = 'r';
    po(j).offColor = [];
    po(j).alphaScale = 0.7;
    po(j).reverse = false;
    po(j).shadedErrorRegion = false;
end
if (nargin == 0)
    hset = po;
    return;
end
if (nargin >= 2 && isstruct(plot_options))
    fn = fieldnames(plot_options);
    
    for j = 1:length(fn)
        for k = 1:(min(length(po), length(plot_options)))
            po(k).(fn{j}) = plot_options(k).(fn{j});
            
        end
    end
end
for j = 1:length(po)
    if (~isempty(po(j).onColor) && ischar(po(j).onColor))
        
        po(j).onColor = char2rgb(po(j).onColor);
    end
    if (~isempty(po(j).offColor) && ischar(po(j).offColor))
        po(j).offColor = char2rgb(po(j).offColor);
    end
    if (~isempty(po(j).backgroundColor) && ischar(po(j).backgroundColor))
        po(j).backgroundColor = char2rgb(po(j).backgroundColor);
        if (j == 1)
            set(0,'DefaultFigureColor', po(j).backgroundColor);
            set(0,'DefaultTextColor', 1 - po(j).backgroundColor);
            set(0,'DefaultAxesColor',  1 - po(j).backgroundColor);
            set(0,'DefaultAxesXColor', 1 - po(j).backgroundColor);
            set(0,'DefaultAxesYColor', 1 - po(j).backgroundColor);
            
        end
    end
end
%po(1)

saveName =  'reorate_vs_time';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))

    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    if (~isempty(po(1).backgroundColor)), set (gcf, 'Color', po(1).backgroundColor); end;
    
    [axunderlay, axplots] = makeUnderlayGraph(ad, po, 'reo_vs_ton', 'reo_vs_toff', 1, po(1).reverse);
    
    if(showtitle),title (axplots, 'Reorientation rate vs. time in cycle');end
    ylabel (axplots, 'rate (min$^{-1}$)');
    xlabel (axplots, 'time (s)');


    emsmallen(axplots, 'FontSize', fontsize, 'Font', font);
end

saveName =  'speed_vs_time';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))

    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    if (~isempty(po(1).backgroundColor)), set (gcf, 'Color', po(1).backgroundColor); end;
    [axunderlay, axplots] = makeUnderlayGraph(ad, po, 'sp_vs_ton', 'sp_vs_toff',60, po(1).reverse);
    if(showlegend), legend (axplots, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation); end
    emsmallen(axplots, 'FontSize', fontsize, 'Font', font);

    if(showtitle),title (axplots, 'Speed vs. time in cycle'); end
    ylabel (axplots, 'speed (cm/m)');
    xlabel (axplots, 'time (s)');
end


saveName = 'reomag_vs_time';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))

    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    if (~isempty(po(1).backgroundColor)), set (gcf, 'Color', po(1).backgroundColor); end;
    [axunderlay, axplots] = makeUnderlayGraph(ad, po, 'reo_mag_vs_ton', 'reo_mag_vs_toff',(180/pi)^2, po(1).reverse);
    if(showlegend), legend (axplots, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation); end

    
%    emsmallen(axplots, 'FontSize', fontsize, 'Font', font);

 %   if(showtitle),title (axplots, 'Size of reorientation vs. time in cycle'); end
 %   ylabel (axplots, 'mean reorientation size (deg)');
    
%     xlabel ('previous heading (degrees)');
    xlabel (axplots, 'time (s)');
    hy = ylabel ('$<\Delta\theta^2>$');
    title ('Mean-Square of reorientations vs. previous heading');
  %  if(showlegend), legend (po.legendEntry, legendOptions{:}); end

    emsmallen(axplots, 'FontSize', fontsize, 'FontName', font);
    yt = get(axplots, 'YTick');
    yt = unique((5*round((sqrt(yt)/5))));
    pos = get(hy, 'Position');
    set (axplots, 'YTick', []);
    xl = get(axplots, 'XLim');
    
    for j = 1:length(yt)
        ht(j) = text(xl(1), yt(j).^2, ['$(' num2str(yt(j)) '^\circ)^2$']);
    end
    set(ht, 'Parent', axplots, 'FontSize', fontsize, 'FontName', font, 'HorizontalAlignment', 'Right');
    set(hy, 'Position', pos);

end

saveName = 'numhs_vs_time';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))

    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    if (~isempty(po(1).backgroundColor)), set (gcf, 'Color', po(1).backgroundColor); end;
    [axunderlay, axplots] = makeUnderlayGraph(ad, po, 'nhs_vs_ton', 'nhs_vs_toff', 1, po(1).reverse);
    if(showlegend), legend (axplots, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation); end
    emsmallen(axplots, 'FontSize', fontsize, 'Font', font);

    if(showtitle),title (axplots, 'number of headsweeps/reorientiation vs. time in cycle'); end
    ylabel (axplots, 'mean num hs/reorientation');
    xlabel (axplots, 'time (s)');
end

saveName = 'val_vs_time';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))

    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    if (~isempty(po(1).backgroundColor)), set (gcf, 'Color', po(1).backgroundColor); end;
    %{
    ax = gca;
    axnew = cloneaxes(ax);

    for j = 1:length(ad)
        newx = (-ad(j).period/2):ad(j).period/2;
        newy = newx;
        inds = 1:ceil((ad(j).period+1)/2);
        inds2 = (length(newy)-length(inds) + 1):length(newy);
        if (po(j).reverse)
            newy(inds) = ad(j).val_vs_toff_highres(inds);
            newy(inds2) = ad(j).val_vs_ton_highres(inds);
        else
            newy(inds2) = ad(j).val_vs_toff_highres(inds);
            newy(inds) = ad(j).val_vs_ton_highres(inds);
        end
        plot (axnew, newx, newy, 'Color', po(j).color, 'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold (axnew,'on');
    end
    
    if(showlegend)
        leghand = legend (axnew, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation);
        if (~isempty(po(1).backgroundColor))
            set(leghand, 'Color', po(1).backgroundColor);
        end
    end
    
    set(axnew, 'XLim', [min(min([ad.etxs]),-max([ad.period])/2), max(max([ad.etxs]), max([ad.period])/2)]);
    props = {'XLim', 'YLim','Position'};
    set(ax, props, get(axnew, props));
    set(axnew, 'Color', 'none');
    set(ax, 'Color', 'none');
    if (po(j).reverse)
        backgroundIntensityUnderlay(ax, ad(1).etxs, ad(1).val_vs_toff,  ad(1).val_vs_ton, po, ad(1).period);
    else
        backgroundIntensityUnderlay(ax, ad(1).etxs, ad(1).val_vs_ton,  ad(1).val_vs_toff, po, ad(1).period);
    end
    set(axnew, 'LineWidth', 1.5);
    axunderlay = ax;
    axplots = axnew;
    %}
    
    [axunderlay, axplots] = makeUnderlayGraph(ad, po, 'val_vs_ton', 'val_vs_toff', 1, po(1).reverse, false);
   
    if(showlegend), legend (axplots, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation); end
    emsmallen(axplots, 'FontSize', fontsize, 'Font', font);

    if(showtitle),title (axplots, [ad(1).fieldname ' vs. time in cycle']); end
    ylabel (axplots, ad(1).fieldname);
    xlabel (axplots, 'time (s)');
    if (isfield(ad, 'rampType') && strcmpi(ad(1).rampType, 'exponential'))
        set(axplots, 'YScale', 'log');
    end
end

exts = {'.tiff', '.eps', '.fig', '.jpg', '.ai', '.pdf'}; ftype = {'tiff', 'eps2c', 'fig', 'jpeg', 'ai', 'pdf'};
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
%------------------------------------------------------------------------%
function nextfigure()
fignum = evalin('caller', 'fignum+1;'); assignin('caller', 'fignum', fignum);
plotNumber = evalin('caller', 'plotNumber+1;'); assignin('caller', 'plotNumber', plotNumber);
forprinting = evalin('caller', 'forprinting;');
if(forprinting)
    figureForPrinting(fignum);
else
    figure(fignum); clf(fignum);
end
set(fignum, 'Color', 'w');  

function backgroundIntensityUnderlay (ax, etxs, fvon, fvoff, po, period)
po = po(1);
if (isempty(po.backgroundColor))
    po.backgroundColor = [1 1 1];
end
%etxs = shifted time axis, fvon fvoff are field vs ton and field vs toff
xl = get(ax, 'XLim');

[x,a] = combineForwardAndReverse (etxs, period, fvon, [], fvoff, []);


ok = x >= xl(1) & x <= xl(2);

x = x(ok);
a = a(ok);
a = a([1 1:end end]);

xl = get(ax, 'XLim');
r = 0.99;
xl = r * xl + (1-r)*xl([2 1]);
x = [xl(1) x xl(2)];
aon = (a - min(a))/(max(a)-min(a));


if (~isempty(po.onColor))
    chigh = po.onColor;
else
    chigh = po.backgroundColor;
end
if (~isempty(po.offColor))
    clow = po.offColor;
else
    clow = po.backgroundColor;
end

cmap = interp1([clow; chigh], linspace(1,2,255));

yl = get(ax, 'YLim');
yl = r * yl + (1-r)*yl([2 1]);
if (~isempty(po.offColor) || ~isempty(po.onColor))
    pcolor (ax,x, yl, [aon;aon]); shading(ax, 'interp'); colormap(ax, cmap); 
end
    
set(ax, 'Color', 'none', 'XTick', [], 'YTick', [],'Layer','bottom','box','off','LineWidth',0.01);

function [axunderlay, axplots] = makeUnderlayGraph (ad, po, fieldon, fieldoff, multiplier, reverse, showerrorbar)
ax = gca;
axnew = cloneaxes(ax);
existsAndDefault('multiplier', 1);
existsAndDefault('reverse', false);
showlegend = evalin('caller','showlegend;');
legendLocation = evalin('caller', 'legendLocation;');
existsAndDefault('showerrorbar', true);

if (~isempty(po(1).backgroundColor))
    bc = po(1).backgroundColor;
    if (ischar(bc))
        bc = char2rgb(bc);
    end
else
    bc = get(gca, 'Color');
end


for j = 1:length(ad)
    if (showerrorbar)
        ebfield_off = ad(j).([fieldoff '_eb']);
        ebfield_on = ad(j).([fieldon '_eb']);
    else
        ebfield_off =zeros(size(ad(j).(fieldoff)));
        ebfield_on = ebfield_off;
    end
    if (reverse)
        [newx, newy, neweb] = combineForwardAndReverse(ad(j).etxs, ad(j).period, multiplier*ad(j).(fieldoff),multiplier*ebfield_off, multiplier*ad(j).(fieldon), multiplier*ebfield_on);
    else
        [newx, newy, neweb] = combineForwardAndReverse(ad(j).etxs, ad(j).period, multiplier*ad(j).(fieldon),multiplier*ebfield_on, multiplier*ad(j).(fieldoff), multiplier*ebfield_off);
    end
    if (showerrorbar)
        if (po(1).shadedErrorRegion)
            axes(axnew);
            
            hh(:,j) = shadedErrorPlot (newx, newy, neweb,[], 'k','LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold (axnew, 'on');
            set(hh(1,j), 'FaceColor', 0.75*bc + 0.25*char2rgb(po(j).color));
            set(hh(2,j), 'Color', char2rgb(po(j).color), 'Marker', po(j).marker);
            
        else
            errorbar (axnew, newx, newy, neweb, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold (axnew,'on');    
        end      
    else
        plot (axnew, newx, newy, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold (axnew,'on');
    end
end

if (po(1).shadedErrorRegion && showerrorbar)
    set(axnew, 'Children', [hh(2,:) hh(1,:)]);
end

if(showlegend)
    leghand = legend (axnew, {po.legendEntry}, 'Interpreter','Latex', 'Location', legendLocation); 
    if (~isempty(po(1).backgroundColor))
        set(leghand, 'Color', po(1).backgroundColor);
    end
end

set(axnew, 'XLim', [min(min([ad.etxs]),-max([ad.period])/2), max(max([ad.etxs]), max([ad.period])/2)]);
props = {'XLim', 'YLim','Position'};
set(ax, props, get(axnew, props));
set(axnew, 'Color', 'none');
set(ax, 'Color', 'none');
if (reverse)
    backgroundIntensityUnderlay(ax, ad(1).etxs, ad(1).val_vs_toff,  ad(1).val_vs_ton, po, ad(1).period);
else
    backgroundIntensityUnderlay(ax, ad(1).etxs, ad(1).val_vs_ton,  ad(1).val_vs_toff, po, ad(1).period);
end

set(axnew, 'LineWidth', 1.5);
axunderlay = ax;
axplots = axnew;

function [newx, newy, neweb] = combineForwardAndReverse (etxs, period, fwd, fwd_eb, rev, rev_eb)

existsAndDefault('fwd_eb', []);
existsAndDefault('rev_eb', []);

doeb = ~isempty(fwd_eb) && ~isempty(rev_eb);
fronteb = [];
backeb = [];
middleeb = [];

if any(etxs == 0)
    middlex = 0;
    middley = 0.5 * (fwd(etxs == 0) + rev(etxs == 0));
    if (doeb)
        middleeb = 0.5 * (fwd_eb(etxs == 0) + rev_eb(etxs == 0));  %not really proper, but oh well, errorbars should be the same anyway
    end
else
    middlex = [];
    middley = [];
  
end

frontx = etxs(etxs > 0)-period/2;
fronty = rev(etxs > 0);

%[frontx, I] = sort (frontx);
%fronty = fronty(I);

if (doeb)
    fronteb = rev_eb(etxs > 0);
 %   fronteb = fronteb(I);
end

backx = etxs(etxs > 0);
backy = fwd(etxs > 0);
if (doeb)
    backeb = fwd_eb(etxs > 0);
end

newx = [frontx middlex backx];
newy = [fronty middley backy];
neweb = [fronteb middleeb backeb];