function figurehandles = spatialNavigationMaggotFiguresSingleExpt (eset, spatial_navigation_options, plot_options, varargin)
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
%   whichGraphs = {}; all if empty - choose from {'DirectionHistogram',
%   'RunStartHistogram','ReorientationRateVsHeading',
%   'InstantaneousDeltaThetaVsTheta','SpeedVsDirection',
%   'ReoDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness',
%   'FirstHeadSwingHandedness', 'ReoDirDistribution','RunLengthHistogram'

SaveDirectory = [];
forprinting = false;
showlegend = false;
whichGraphs = {};
startFigNum = 1;
backgroundColor = [];
legendLocation = 'BestOutside';
showtitle = true;
vectorgraphics = [];
if (isempty(vectorgraphics))
    vectorgraphics = ~(isempty(SaveDirectory) && ~showlegend);
end

varargin = assignApplicable(varargin);
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


ccc = 'bgrcymk';
sss = 'sodvh>p^<';
po.lineWidth = 2;
po.color = ccc(1);

po.marker = sss(1);
po.plotOptions = {};
po.shadedErrorRegion = false;
po.useGauss = false;
po.directions = [-180,-90, 0, 90];
po.colors = {[0.1 0 1], [0.1 0.8 0.2], [1 0 0], [0.6 0.6 0.0]};

if (nargin == 0)
    figurehandles = po;
    return;
end

if (nargin > 1 && isstruct(spatial_navigation_options))
    fn = fieldnames(spatial_navigation_options);
    for j = 1:length(fn)
        sno.(fn{j}) = spatial_navigation_options.(fn{j});
    end
end

if (isa(eset, 'ExperimentSet'))
    ad = spatialNavigationMaggotAnalysis(eset, sno);
else
    ad = eset;
end

nexp = length(ad);
if (nexp > 1)
    error ('call on only one experiment analyzed data set');
end

po.legendEntry = cellfun(@(x) ['to ' num2str(x) '$^\circ$'], num2cell(po.directions),'UniformOutput',false);
if (nargin > 2 && isstruct(plot_options))
    fn = fieldnames(plot_options);
    for j = 1:length(fn)
        po.(fn{j}) = plot_options.(fn{j});
    end
end
for j = 1:length(po.colors)
    if (ischar(po.colors{j}))
        po.colors{j} = char2rgb(po.colors{j});
    end
end


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
    
   
    quadrantGraph(ad, 'txc', 'thetahist', po, 1/sum(ad.thetahist), true);
%     
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).thetahist([1:end 1]), ad(j).thetahist_eb([1:end 1]), 'Color', po(j).color, 'Marker', po(j).marker, 'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end

    xlabel ('instantaneous run heading (degrees)');
    ylabel ('relative frequency');
    title ('Instantaneous Direction of Forward Movement');
    if (showlegend), legend (po.legendEntry, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end

end

saveName =  'ReorientationRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    quadrantGraph(ad, 'txc', 'reohist', po, 1, true);
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).reohist([1:end 1]),  ad(j).reohist_eb([1:end 1]), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('reorientation rate (min$^{-1}$)');
    title ('Reorientation Rate vs. Direction of Forward Movement');
    if(showlegend), legend (po.legendEntry, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'RunStartHistogram';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    quadrantGraph(ad, 'txc', 'runStartDirectionHist', po, 1, true);
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).runStartDirectionHist([1:end 1]),ad(j).runStartDirectionHist_eb([1:end 1]),'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('start direction of run (degrees)');
    ylabel ('relative frequency');
    title ('Heading at beginning of runs');
    if(showlegend), legend (po.legendEntry,legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName = 'RunLengthHistogram';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    
    maxI = find(cumsum(ad.runTimeHistTowards) > 0.90, 1, 'first');
    maxI = max(maxI, find(cumsum(ad.runTimeHistAway) > 0.90, 1, 'first'));
    inds = 1:maxI;
    
    [~,I] = max(cosd(po.directions - sno.preferredDirection));
    c = po.colors{I};
    errorbar(ad.runTimeAxis(inds), ad.runTimeHistTowards(inds),ad.runTimeHistTowards_eb(inds), 'ko-', 'Color', c, 'LineWidth', 2); hold on
    mrl = mean(ad.runTime(cos(ad.runMeanTheta - deg2rad(sno.preferredDirection)) > 1/sqrt(2) & ad.runTime < ad.runTimeAxis(maxI)));
   % plot (ad.runTimeAxis(inds), diff(ad.runTimeAxis(1:2))/mrl.*exp(-ad.runTimeAxis(inds)./mrl), 'k--', 'Color', c);
    [~,I2] = min(cosd(po.directions - sno.preferredDirection));
    c = po.colors{I2};
    errorbar(ad.runTimeAxis(inds), ad.runTimeHistAway(inds),ad.runTimeHistAway_eb(inds), 'ks-', 'Color', c, 'LineWidth', 2); hold on
    mrl = mean(ad.runTime(-cos(ad.runMeanTheta - deg2rad(sno.preferredDirection)) > 1/sqrt(2) & ad.runTime < ad.runTimeAxis(maxI)));
    %plot (ad.runTimeAxis(inds), diff(ad.runTimeAxis(1:2))/mrl.*exp(-ad.runTimeAxis(inds)./mrl), 'k--', 'Color', c);
    set(gca, 'YScale', 'log');
    
    
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).runStartDirectionHist([1:end 1]),ad(j).runStartDirectionHist_eb([1:end 1]),'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('run duration (s)');
    ylabel ('relative frequency');
    title ('histogram of run lengths');
    le = cellfun(@(x) ['run heading towards ' num2str(x) '$^\circ$'], num2cell(po.directions([I I2])),'UniformOutput',false);
    if(showlegend), legend (le,legendOptions{:}); end
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
    if (po.useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'txc';
        iscirc = true;
    end
    
    quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
%     for j = 1:nexp
%         errorbar (ad(j).txc, rad2deg(ad(j).instantaneousdthetavstheta([1:end 1])), rad2deg(ad(j).instantaneousdthetavstheta_eb([1:end 1])),'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    title ('Heading Change vs. Direction (in runs)');
    xlabel ('heading (degrees)');
    ylabel ('mean heading change (degrees/min)');
    if(showlegend), legend (po.legendEntry, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'SpeedVsDirection';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    yfield = 'speedVsDir';
    ymult = 60;
    if (po.useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'txc';
        iscirc = true;
    end
    
    quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
    %quadrantGraph(ad, 'speedVsDir', po, 60);
%     for j = 1:nexp
%         errorbar (ad(j).txc, ad(j).speedVsDir([1:end 1]), ad(j).speedVsDir_eb([1:end 1]), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('heading (degrees)');
    ylabel ('speed (cm/min)');
    title ('Speed vs. Direction (in runs)');
    if(showlegend), legend (po.legendEntry, legendOptions{:}); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end
dtxf = -180:10:180;

[rdfit, rmfit, ~, rdfitci, rmfitci] = ad.reo_dir_model.reorientationMeansVsTheta(ad.reo_dir_model.params, deg2rad(dtxf), ad.reo_dir_model.hessian, 0.95);

rdfit = rad2deg(rdfit);
rdfitci = rad2deg(rdfitci);
rmfit = rmfit*(180/pi)^2;
rmfitci = rmfitci*(180/pi)^2;
saveName =  'ReoDirVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    yfield = 'reoDir';
    ymult = 180/pi;
    if (po.useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'reotxc';
        iscirc = true;
    end
    
    quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
%    quadrantGraphReo(ad, 'reoDir', po, 180/pi);
    hold on;
    plot (dtxf, rdfit, 'k--', dtxf, rdfitci(1,:), 'k:', dtxf, rdfitci(2,:), 'k:', 'LineWidth', 2);
    hold off;
    
%     for j = 1:nexp
%         errorbar (ad(j).reotxc, rad2deg(ad(j).reoDir([1:end 1])), rad2deg(ad(j).reoDir_eb([1:end 1])), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('previous heading (degrees)');
    ylabel ('mean reorientation direction');
    title ('Mean reorientations vs. previous heading');
    if(showlegend), legend (po.legendEntry, legendOptions{:}); end

    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
end

saveName =  'ReoDirDistribution';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
   
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    %po.directions = [-180,-90, 0, 90];
    %po.colors = {[0 0 1], [0 1 0], [1 0 0], [1 1 0]};
    for j = 1:4
        
        spx(j) = subplot(2,2,j);
        [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(j) + 180, 360)- 180));
        denom = sum(ad.reo_dtheta_dist_hs{j});% + ad.reo_dtheta_dist_nohs{j});
        if (isempty(backgroundColor))
            bc = [1 1 1];
        else
            bc = backgroundColor;
        end
        fc1 = bc;
        fc2 = 0.4*po.colors{ind} + 0.6*bc;
        
        bardat = [ad.reo_dtheta_dist_hs{j}([1:end 1])]/denom;%;ad.reo_dtheta_dist_nohs{j}([1:end 1])]'/denom;
        hh = bar(ad.dtxc,  bardat, 'FaceColor', fc2, 'EdgeColor', po.colors{ind}, 'LineWidth', 1); hold on
        %set(hh(1), 'FaceColor', fc2);
%        errorbar(ad.dtxc, ad.reo_dtheta_dist{j}([1:end 1])/denom,ad.reo_dtheta_dist_eb{j}([1:end 1])/denom, 'k.', 'LineWidth', 2, 'Color', po.colors{ind});
        plot (ad.dtxf, ad.reo_dtheta_fit{j}/denom, 'k-','LineWidth',3,'Color', po.colors{ind});
        plot (ad.dtxf, ad.reo_dtheta_fit_r{j}/denom, 'k--', ad.dtxf, ad.reo_dtheta_fit_l{j}/denom, 'k--', 'LineWidth',2,'Color', po.colors{ind});
        xlabel ('change in heading (deg)');
        ylabel ('fraction of reorientations');
        title (['previous heading near ' num2str(ad.reobasedirections(j)) '$^\circ$']);
        emsmallen(gca, 'FontSize', fontsize, 'Font', font);
        if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
        hold off
        
    end
    %{
    hold (spx(1), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(2) + 180, 360)- 180));
    plot (spx(1), ad.dtxf, ad.reo_dtheta_fit{2}/sum(ad.reo_dtheta_dist{2}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(1), 'off');
    hold (spx(2), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(1) + 180, 360)- 180));
    plot (spx(2), ad.dtxf, ad.reo_dtheta_fit{1}/sum(ad.reo_dtheta_dist{1}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(2), 'off');
    %}
    
    
    disp('3');
    
    ym = 0;
    for j = 1:4
       ym = max(ym, max(get(spx(j), 'YLim')));
    end
    [~,I] = max(ad.reo_dtheta_fit{1}.*(ad.dtxf < 0));
    [~,J] = max(ad.reo_dtheta_fit{1}.*(ad.dtxf > 0));
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(1) + 180, 360)- 180));
    towardscol = po.colors{ind};
    towardspeaks = ad.dtxf([I J]);
    disp('4')
    %{
    [~,I] = max(ad.reo_dtheta_fit{2}.*(ad.dtxf < 0));
    [~,J] = max(ad.reo_dtheta_fit{2}.*(ad.dtxf > 0));
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(2) + 180, 360)- 180));
    awaycol = po.colors{ind};
    awaypeaks = ad.dtxf([I J]);
    hold (spx(1), 'on');
    plot (spx(1), towardspeaks(1)*[1 1], [0 ym], 'k--', towardspeaks(2)*[1 1], [0 ym], 'k--', 'Color', towardscol,'LineWidth',2);
    plot (spx(1), awaypeaks(1)*[1 1], [0 ym], 'k--', awaypeaks(2)*[1 1], [0 ym], 'k--', 'Color', awaycol,'LineWidth',2);
    hold (spx(1), 'off');
    hold (spx(2), 'on');
    plot (spx(2), towardspeaks(1)*[1 1], [0 ym], 'k--', towardspeaks(2)*[1 1], [0 ym], 'k--', 'Color', towardscol,'LineWidth',2);
    plot (spx(2), awaypeaks(1)*[1 1], [0 ym], 'k--', awaypeaks(2)*[1 1], [0 ym], 'k--', 'Color', awaycol,'LineWidth',2);
    hold (spx(2), 'off');
    %}
    set(spx, 'YLim', [0 ym]);
    
end


saveName =  'ReoDirDistributionPolar';
if (false && isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    for j = 1:4
        spx(j) = subplot(2,2,j);
        [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(j) + 180, 360)- 180));
        denom = sum(ad.reo_dtheta_dist_hs{j} + ad.reo_dtheta_dist_nohs{j});
        if (isempty(backgroundColor))
            bc = [1 1 1];
        else
            bc = backgroundColor;
        end
            fc1 = bc;
        	fc2 = 0.4*po.colors{ind} + 0.6*bc;
        bardat = [ad.reo_dtheta_dist_nohs{j}([1:end 1]); ad.reo_dtheta_dist_hs{j}([1:end 1])]'/denom;
        hh = bar(ad.dtxc,  bardat, 'Stacked', 'FaceColor', fc2, 'EdgeColor', po.colors{ind}, 'LineWidth', 1); hold on
        set(hh(1), 'FaceColor', fc1);
%        errorbar(ad.dtxc, ad.reo_dtheta_dist{j}([1:end 1])/denom,ad.reo_dtheta_dist_eb{j}([1:end 1])/denom, 'k.', 'LineWidth', 2, 'Color', po.colors{ind});
        plot (ad.dtxf, ad.reo_dtheta_fit{j}/denom, 'k-','LineWidth',3,'Color', po.colors{ind});
        plot (ad.dtxf, ad.reo_dtheta_fit_r{j}/denom, 'k--', ad.dtxf, ad.reo_dtheta_fit_l{j}/denom, 'k--', 'LineWidth',2,'Color', po.colors{ind});
        xlabel ('change in heading (deg)');
        ylabel ('fraction of reorientations');
        title (['previous heading near ' num2str(ad.reobasedirections(j)) '$^\circ$']);
        emsmallen(gca, 'FontSize', fontsize, 'Font', font);
        if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
        hold off
    end
    %{
    hold (spx(1), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(2) + 180, 360)- 180));
    plot (spx(1), ad.dtxf, ad.reo_dtheta_fit{2}/sum(ad.reo_dtheta_dist{2}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(1), 'off');
    hold (spx(2), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(1) + 180, 360)- 180));
    plot (spx(2), ad.dtxf, ad.reo_dtheta_fit{1}/sum(ad.reo_dtheta_dist{1}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(2), 'off');
    %}
    
    
    
    
    ym = 0;
    for j = 1:4
       ym = max(ym, max(get(spx(j), 'YLim')));
    end
    [~,I] = max(ad.reo_dtheta_fit{1}.*(ad.dtxf < 0));
    [~,J] = max(ad.reo_dtheta_fit{1}.*(ad.dtxf > 0));
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(1) + 180, 360)- 180));
    towardscol = po.colors{ind};
    towardspeaks = ad.dtxf([I J]);
    %{
    [~,I] = max(ad.reo_dtheta_fit{2}.*(ad.dtxf < 0));
    [~,J] = max(ad.reo_dtheta_fit{2}.*(ad.dtxf > 0));
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(2) + 180, 360)- 180));
    awaycol = po.colors{ind};
    awaypeaks = ad.dtxf([I J]);
    hold (spx(1), 'on');
    plot (spx(1), towardspeaks(1)*[1 1], [0 ym], 'k--', towardspeaks(2)*[1 1], [0 ym], 'k--', 'Color', towardscol,'LineWidth',2);
    plot (spx(1), awaypeaks(1)*[1 1], [0 ym], 'k--', awaypeaks(2)*[1 1], [0 ym], 'k--', 'Color', awaycol,'LineWidth',2);
    hold (spx(1), 'off');
    hold (spx(2), 'on');
    plot (spx(2), towardspeaks(1)*[1 1], [0 ym], 'k--', towardspeaks(2)*[1 1], [0 ym], 'k--', 'Color', towardscol,'LineWidth',2);
    plot (spx(2), awaypeaks(1)*[1 1], [0 ym], 'k--', awaypeaks(2)*[1 1], [0 ym], 'k--', 'Color', awaycol,'LineWidth',2);
    hold (spx(2), 'off');
    %}
    set(spx, 'YLim', [0 ym]);
    
end

saveName =  'HsDirDistribution';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    for j = 1:4
        spx(j) = subplot(2,2,j);
        [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(j) + 180, 360)- 180));
        denom = sum(ad.hs_all_dist{j});
        if (~isempty(backgroundColor))
           fc = 0.2 * po.colors{ind} + 0.8 * backgroundColor;
           fc2 = backgroundColor;
        else
           fc = 0.2*po.colors{ind} + 0.8*[1 1 1];
           fc2 = [1 1 1];
        end
        hh = bar(ad.hsdtx, [ad.hs_acc_dist{j};ad.hs_rej_dist{j}]'/denom, 'Stacked');
        set (hh(1), 'FaceColor', fc, 'EdgeColor', po.colors{ind}, 'LineWidth', 1); 
        set (hh(2), 'FaceColor', fc2, 'EdgeColor', po.colors{ind}, 'LineWidth', 1);
        xlabel ('body bend angle (deg)');
        ylabel ('fraction of headswings');
        title (['previous heading near ' num2str(ad.reobasedirections(j)) '$^\circ$']);
%        legend ('accepted', 'rejected', 'Location', 'BestOutside');
        emsmallen(gca, 'FontSize', fontsize, 'Font', font);
        if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
        hold off
    end
    %{
    hold (spx(1), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(2) + 180, 360)- 180));
    plot (spx(1), ad.dtxf, ad.reo_dtheta_fit{2}/sum(ad.reo_dtheta_dist{2}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(1), 'off');
    hold (spx(2), 'on');
    [~,ind] = min(abs(mod(po.directions - ad.reobasedirections(1) + 180, 360)- 180));
    plot (spx(2), ad.dtxf, ad.reo_dtheta_fit{1}/sum(ad.reo_dtheta_dist{1}), 'k--','LineWidth',2,'Color', po.colors{ind});
    hold (spx(2), 'off');
    %}
    
    
    
    
    ym = 0;
    for j = 1:4
       ym = max(ym, max(get(spx(j), 'YLim')));
    end
   
    set(spx, 'YLim', [0 ym]);
    
end

saveName =  'ReoMagVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

     yfield = 'reoMag';
    ymult = (180/pi)^2;
    if (po.useGauss)
        xfield = 'txf';
        yfield = [yfield '_gauss'];
        iscirc = false;
    else
        xfield = 'reotxc';
        iscirc = true;
    end
    
    quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
    hold on;
    plot (dtxf, rmfit, 'k--', dtxf, rmfitci(1,:),'k:', dtxf, rmfitci(2,:), 'k:','LineWidth', 2);
    hold off;
    
%     for j = 1:nexp
%         errorbar (ad(j).reotxc, rad2deg(ad(j).reoMag([1:end 1])), rad2deg(ad(j).reoStd([1:end 1])), 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
%     end
    xlabel ('previous heading (degrees)');
    hy = ylabel ('$<\Delta\theta^2>$');
    title ('Mean-Square of reorientations vs. previous heading');
    if(showlegend), legend (po.legendEntry, legendOptions{:}); end

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
    if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); set(ht, 'Color', 1-backgroundColor); end
end

saveName =  'HeadSwingAcceptanceHandedness';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    
    hsbar = zeros(1, length(ad.hstx));
    hseb =  zeros(1, length(ad.hstx));
    
    hsbar (1, :) = ad.headSwingAcceptanceRateLeft;
    hsbar (2, :) = ad.headSwingAcceptanceRateRight;
    hseb (1, :) = ad.headSwingAcceptanceRateLeft_eb; 
    hseb (2, :) = ad.headSwingAcceptanceRateRight_eb;
    
    
    gnames = cellfun(@num2str, num2cell(ad.hstx),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    ttl = ('Head Swing Acceptances By Bend Direction');
    
    le = {'to left', 'to right'};
   
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl,[],[],le,[],'axis'); 
    for j = 1:4
        I = find(cosd(ad.hstx(j)) == cosd(po.directions) & sind(ad.hstx(j)) == sind(po.directions), 1);
        if (~isempty(I))
            set(hhh.bars(j,1), 'FaceColor', po.colors{I}, 'EdgeColor', po.colors{I});
            set(hhh.bars(j,2), 'FaceColor', 'w', 'EdgeColor',po.colors{I}, 'LineWidth', 4);
        end
    end
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    if (showtitle)
        title (ttl);
          p = get(gca, 'Position');
        p(4) = 0.95 * p(4); set(gca, 'Position', p);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);

    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end

%{
saveName =  'HeadSwingAcceptanceDirection';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;


   
    hsbar = zeros(2*length(ad), length(ad.hstx));
    hseb =  zeros(2*nexp, length(ad.hstx));
    
    hsbar (1, :) = ad.headSwingAcceptanceRateTowards;
    hsbar (2, :) = ad.headSwingAcceptanceRateAway;
    hseb (1, :) = ad.headSwingAcceptanceRateTowards_eb; 
    hseb (2, :) = ad.headSwingAcceptanceRateAway_eb;
    gnames = cellfun(@num2str, num2cell(ad.hstx),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
     le = {['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['away from ' num2str(sno.preferredDirection) '$^{\circ}$']};
    ttl = ('Head Swing Acceptances By Bend Direction');
    
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl,[],[],le,[],'axis'); 
    for j = 1:4
        I = find(cosd(ad.hstx(j)) == cosd(po.directions) & sind(ad.hstx(j)) == sind(po.directions), 1);
        if (~isempty(I))
            set(hhh.bars(j,1), 'FaceColor', po.colors{I}, 'EdgeColor', po.colors{I});
            set(hhh.bars(j,2), 'FaceColor', 'w', 'EdgeColor',po.colors{I}, 'LineWidth', 4);
        end
    end
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
    if (showtitle)
        title (ttl);
          p = get(gca, 'Position');
        p(4) = 0.95 * p(4); set(gca, 'Position', p);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    
    
%     hsbar = zeros(nexp, 2);
%     hsbar(:, 1) = [ad.headSwingAcceptanceRateTowards];
%     hsbar(:, 2) = [ad.headSwingAcceptanceRateAway];
%     hseb =  hsbar;
%     hseb(:, 1) = [ad.headSwingAcceptanceRateTowards_eb];
%     hseb(:, 2) = [ad.headSwingAcceptanceRateAway_eb];
%     gnames = po.legendEntry;
%     ylbl = 'probability of accepting head swing';
%     ttl = ('Head Swing Acceptances By Bend Direction');
%     bardat = hsbar;
%     ebdat = hseb;
%     hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],[],[]);
%     th = rotateticklabel(hhh.ax, 20);
%     set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
%     set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['away from ' num2str(sno.preferredDirection) '$^{\circ}$']),legendOptions{:});
%     set(hhh.ax, 'YLim', [0 1]);
%      pos = get(hhh.ax, 'Position');
%     w = pos(3);
%     h = pos(4);
%     r = 0.15;
%     pos = pos + [w*r/2 h*r/2 -w*r -h*r];
%     set(hhh.ax, 'Position', pos);
%     if (isempty(backgroundColor))
%         set (hhh.bars(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', [0.7 0.7 0.7]);
%         set (hhh.bars(2), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 3);
%     else
%         set (hhh.bars(1), 'FaceColor', 0.7 - 0.7*backgroundColor, 'EdgeColor', 0.3 - 0.3*backgroundColor);
%         set (hhh.bars(2), 'FaceColor', 0.3 - 0.3*backgroundColor, 'EdgeColor', 0.7 - 0.7*backgroundColor, 'LineWidth', 3);
%     end
%     emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
     if (~isempty(backgroundColor)), set (gca, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
     if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end
%}

%{
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
%}
saveName =  'FirstHeadSwingHandedness';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    gnames = cellfun(@num2str, num2cell(ad(1).hstx),'UniformOutput',false);
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;
    ttl = 'First Head Swing Bias vs. Previous Direction';
    ylbl = 'probability first head sweep is to left';
    le = {po(1:nexp).legendEntry};
    
    
    gnames = cellfun(@num2str, num2cell(ad.hstx),'UniformOutput',false);
    ylbl = 'probability first head sweep is to left';
  
    ttl = ('Head Swing Initiation Bias');
    bardat = ad.firstHSDir';
    ebdat = ad.firstHSDir_eb';

    hhh = barweb_marc(bardat, ebdat, 0.8, gnames,[],  [],ylbl);
    for j = 1:4
        I = find(cosd(ad.hstx(j)) == cosd(po.directions) & sind(ad.hstx(j)) == sind(po.directions), 1);
        if (~isempty(I))
            set(hhh.bars(j), 'FaceColor', po.colors{I}, 'EdgeColor', po.colors{I});
        end
    end
    if (showtitle)
        title (ttl);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    hold on;
    plot (get (gca, 'XLim'), [0.5 0.5], 'r--', 'LineWidth', 3);
    set(hhh.baseline, 'LineWidth', 4);
    
%     %{
%     bardat = firstHSDir';
%     ebdat = firstHSDir_eb';
%     %}
%   
%     %Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
%     if (nexp > 1)
%         hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
%     else
%         hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
%       %  set(hhh.ax, 'XTick', 1:length(bardat), 'XTickLabel', gnames);
%     end
%     for j = 1:nexp
%         set(hhh.bars(j), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
%     end
%     %xlabel ('previous direction');
%     %ylabel ('probability first head sweep is to left');
%     %title ('First Head Sweep Directional Bias');
%     emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
%     if (nexp > 1)
%         set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
%     end
%     hold on
%     xl = get(gca, 'Xlim');
%     plot (xl, [0.5 0.5], 'r--','LineWidth', 3);
%     hold off
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end


%{
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
    gnames = po.legendEntry;
    ylbl = '$\bar{v}/\bar{s}$';
    ttl = ('Navigation Indices');
    bardat = hsbar;
    ebdat = hseb;
    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl);
    th = rotateticklabel(hhh.ax, 20);
    set(th, 'Interpreter', 'Latex', 'FontSize', fontsize-2, 'FontName', font);
    
    set(legend(['toward ' num2str(sno.preferredDirection) '$^{\circ}$'], ['toward ' num2str(mod(sno.preferredDirection + 90, 360)) '$^{\circ}$']), legendOptions{:});
 %   set(hhh.ax, 'YLim', [-1 1]);
    
    if (isempty(backgroundColor))
        set (hhh.bars(1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', [0.7 0.7 0.7]);
        set (hhh.bars(2), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 3);
    else
        set (hhh.bars(1), 'FaceColor', 0.7 - 0.7*backgroundColor, 'EdgeColor', 0.3 - 0.3*backgroundColor);
        set (hhh.bars(2), 'FaceColor', 0.3 - 0.3*backgroundColor, 'EdgeColor', 0.7 - 0.7*backgroundColor, 'LineWidth', 3);
    end
    emsmallen(gca, 'FontSize', fontsize, 'FontName', font);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end



saveName =  'StrategicIndices';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    

    hsbar = [[ad.reoRateIndex];[ad.reoMagIndex];[ad.reoDirIndex];[ad.headSwingAcceptanceRateIndex]];
    hseb = [[ad.reoRateIndex_eb];[ad.reoMagIndex_eb];[ad.reoDirIndex_eb];[ad.headSwingAcceptanceRateIndex_eb]];
    gnames = {{'reorientation','rate'}, {'reorientation','magnitude'},{'reorientation','direction'} {'head swing', 'acceptance'}};
    le = po.legendEntry;
    ylbl = 'index';
    ttl = ('Navigation Indices');
    bardat = hsbar;
    ebdat = hseb;
    hhh = barweb(bardat, ebdat, 0.8, [],ttl,  [],ylbl,[],[],le,[],'axis');
    
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
    pos = get(hhh.ax, 'Position');
    w = pos(3);
    h = pos(4);
    r = 0.15;
    pos = pos + [w*r/2 h*r/2 -w*r -h*r];
    set(hhh.ax, 'Position', pos);
    if (~isempty(backgroundColor)), set (hhh.ax, 'Color', backgroundColor, 'XColor', 1- backgroundColor, 'YColor', 1 - backgroundColor); end
    if (~isempty(backgroundColor)), set (hhh.errors, 'Color', 1-backgroundColor); end
end
%}

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


figurehandles = hset;
set(0, default_properties);

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

    
function quadrantGraph(ad, xfield, yfield,po, mult, iscirc)
%function quadrantGraph(ad, xfield, yfield,po, mult, iscirc)
[directions,I] = sort(po.directions);
colors = po.colors(I);
existsAndDefault('mult', 1);

tx = ad.(xfield);


q = zeros(size(tx));
for j = 1:length(tx)
    [~,I] = max(cosd(directions - tx(j)));
    q(j) = I;
end
yd = mult*ad.(yfield);
ebd = mult*ad.([yfield '_eb']);

if (iscirc)
    yd = yd([1:end 1]);
    ebd = ebd([1:end 1]);
end

for j = 1:length(directions)
    inds = find(q == j);
    [xdata{j},I] = sort(tx(inds));
    inds = inds(I);
    ydata{j} = yd(inds);
    ebdata{j} = ebd(inds);
end
for j = 1:length(xdata)
    if (any(diff(xdata{j}) > min(diff(directions))))        
        [~,I] = max(diff(xdata{j}));
        ni = length(xdata) + 1;
        colors{ni} = colors{j};
        xdata{ni} = xdata{j}((I+1):end);
        xdata{j} = xdata{j}(1:I);        
        ydata{ni} = ydata{j}((I+1):end);
        ydata{j} = ydata{j}(1:I);
        ebdata{ni} = ebdata{j}((I+1):end);
        ebdata{j} = ebdata{j}(1:I);
    end
end



xf = zeros(size(xdata)); xr = xf; yf = xf; yr = xf; ebf = xf; ebr = xf;

for j = 1:length(xdata)
    ahead = mod(j, length(xdata)) + 1;
    behind = mod(j-2, length(xdata)) + 1;
   
    xf(j) = 0.5* (xdata{j}(end) + mod(xdata{ahead}(1) - xdata{j}(end), 360) + xdata{j}(end));
    xr(j) = 0.5* (-mod(xdata{j}(1) - xdata{behind}(end), 360) + xdata{j}(1) + xdata{j}(1));
    yf(j) = 0.5* (ydata{j}(end) + ydata{ahead}(1));
    yr(j) = 0.5* (ydata{behind}(end) + ydata{j}(1));
    ebf(j) = 0.5* (ebdata{j}(end) + ebdata{ahead}(1));
    ebr(j) = 0.5* (ebdata{behind}(end) + ebdata{j}(1));
end

if (po.shadedErrorRegion)
    for j = 1:length(xdata)
        xxdata{j} = [xr(j) xdata{j} xf(j)];
        yydata{j} = [yr(j) ydata{j} yf(j)];
        eebdata{j} = [ebr(j) ebdata{j} ebf(j)];
    end
    hh = shadedErrorPlot(xxdata, yydata, eebdata, [], colors, 'LineWidth', po.lineWidth);
    set(hh((length(xdata)+1):end), 'LineWidth', po.lineWidth, po.plotOptions{:});
    hold on;
    if (length(tx) > 60)
        marker = 'none';
    else
        marker = po.marker;
    end
    for j = 1:length(xdata)
        plot (xdata{j}, ydata{j}, 'Color', colors{j},'Marker', marker, 'LineStyle', 'none', 'LineWidth', po.lineWidth, po.plotOptions{:}); 
    end
    hold off;
else
    for j = 1:length(xdata)
        plot ([xr(j) xdata{j} xf(j)], [yr(j) ydata{j} yf(j)], 'k-', 'Color', colors{j}, 'LineWidth', po.lineWidth, po.plotOptions{:}); hold on;
        errorbar (xdata{j}, ydata{j}, ebdata{j}, 'Color', colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth, po.plotOptions{:});
    end
    hold off
end

%     
% function quadrantGraph(ad, fieldname,po, mult)
% %function quadrantGraph(ad, fieldname,po, mult)
% directions = po.directions;
% colors = po.colors;
% existsAndDefault('mult', 1);
% q = zeros(size(ad.txc));
% for j = 1:length(ad.txc)
%     [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
%     q(j) = I;
% end
% 
% yd = mult*ad.(fieldname)([1:end 1]);
% ebd = mult*ad.([fieldname '_eb'])([1:end 1]);
% 
% 
% for j = 1:(length(ad.txc) - 1)
%     xdata = interp1(ad.txc, [j j+0.5]);
%     ydata = interp1(yd, [j j+0.5]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
%     xdata = interp1(ad.txc, [j+0.5 j+1]);
%     ydata = interp1(yd, [j+0.5 j+1]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end
% 
% for j = 1:length(directions)
%     inds = q == j;
%     
%     hhhh(j) = errorbar (ad.txc(inds), yd(inds), ebd(inds), 'k.', 'Color',colors{j}); hold on;
%    % get(hhhh)
%    % bob = get(get(hhhh,'Children'))
% %    pause
% end
% set(hhhh, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:});
% 
% function quadrantGraphReo(ad, fieldname,po, mult)
% %function quadrantGraph(ad, fieldname,po, mult)
% directions = po.directions;
% colors = po.colors;
% existsAndDefault('mult', 1);
% q = zeros(size(ad.reotxc));
% for j = 1:length(ad.reotxc)
%     [~,I] = max(cosd(directions).*cosd(ad.reotxc(j)) + sind(directions).*sind(ad.reotxc(j)));
%     q(j) = I;
% end
% 
% yd = mult*ad.(fieldname)([1:end 1]);
% ebd = mult*ad.([fieldname '_eb'])([1:end 1]);
% 
% for j = 1:length(directions)
%     inds = q == j;
%     errorbar (ad.reotxc(inds), yd(inds), ebd(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end
% 
% for j = 1:(length(ad.reotxc) - 1)
%     xdata = interp1(ad.reotxc, [j j+0.5]);
%     ydata = interp1(yd, [j j+0.5]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
%     xdata = interp1(ad.reotxc, [j+0.5 j+1]);
%     ydata = interp1(yd, [j+0.5 j+1]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end