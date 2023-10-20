function hset = checkerNavigationFigures(ad, plot_options, varargin)
%function hset = checkerNavigationFigures(ad, plot_options, varargin)
%'HeadSwingAcceptanceHandedness_0', 'HeadSwingAcceptanceDirection_0',
%'HeadSwingAcceptanceHandedness_45','HeadSwingDirection_0,
%'HeadSwingDirection_45', 'ReorientationRateVsHeading',
%'ReorientationRateVsDistance', 'ReorientationMagVsHeading',
%'ReorientationDirVsHeading', 'PauseRateVsHeading'


if (nargin > 0 && length(ad) == 1)
    hset = checkerNavigationFiguresSingleExpt(ad, plot_options, varargin{:});
    return;
end

SaveDirectory = [];
forprinting = false;
showlegend = false;
fignum = 0;
varargin = assignApplicable(varargin);

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

saveName =  'HeadSwingAcceptanceHandedness_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx0));
    hseb =  zeros(2*nexp, length(ad.hstx0));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_acctoleft0; 
        hsbar (2:2:end, :) = ad(j).hs_acctoright0;
        hseb (1:2:end, :) = ad(j).hs_acctoleft0_eb; 
        hseb (2:2:end, :) = ad(j).hs_acctoright0_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
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
end
saveName =  'HeadSwingAcceptanceDirection_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx0));
    hseb =  zeros(2*nexp, length(ad.hstx0));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_acctolight0; 
        hsbar (2:2:end, :) = ad(j).hs_acctodark0;
        hseb (1:2:end, :) = ad(j).hs_acctolight0_eb; 
        hseb (2:2:end, :) = ad(j).hs_acctodark0_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' to light'];
        le{(2*j)} = [po(j).legendEntry ' to dark'];
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
end

saveName =  'HeadSwingAcceptanceHandedness_45';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx45));
    hseb =  zeros(2*nexp, length(ad.hstx45));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_acctoleft45; 
        hsbar (2:2:end, :) = ad(j).hs_acctoright45;
        hseb (1:2:end, :) = ad(j).hs_acctoleft45_eb; 
        hseb (2:2:end, :) = ad(j).hs_acctoright45_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx45),'UniformOutput',false);
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
end

saveName =  'HeadSwingAcceptanceDirection_45';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx45));
    hseb =  zeros(2*nexp, length(ad.hstx45));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_acctolight45; 
        hsbar (2:2:end, :) = ad(j).hs_acctodark45;
        hseb (1:2:end, :) = ad(j).hs_acctolight45_eb; 
        hseb (2:2:end, :) = ad(j).hs_acctodark45_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx45),'UniformOutput',false);
    ylbl = 'probability of accepting head swing';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' to light'];
        le{(2*j)} = [po(j).legendEntry ' to dark'];
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
end

saveName =  'HeadSwingDirection_0';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx0));
    hseb =  zeros(2*nexp, length(ad.hstx0));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_dir0; 
        hsbar (2:2:end, :) = ad(j).firsths_dir0;
        hseb (1:2:end, :) = ad(j).hs_dir0_eb; 
        hseb (2:2:end, :) = ad(j).firsths_dir0_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx0),'UniformOutput',false);
    ylbl = 'probability head sweep is to left';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' all hs'];
        le{(2*j)} = [po(j).legendEntry ' first hs'];
    end
    ttl = ('Head Swing Initiation Bias');
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    for j = 1:nexp
        set(hhh.bars(2*j - 1), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
        set(hhh.bars(2*j), 'FaceColor', 'w', 'EdgeColor', po(j).color, 'LineWidth', 4);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
end

saveName =  'HeadSwingDirection_45';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    hsbar = zeros(2*length(ad), length(ad.hstx45));
    hseb =  zeros(2*nexp, length(ad.hstx45));
    for j = 1:length(ad)
        hsbar (1:2:end, :) = ad(j).hs_dir45; 
        hsbar (2:2:end, :) = ad(j).firsths_dir45;
        hseb (1:2:end, :) = ad(j).hs_dir45_eb; 
        hseb (2:2:end, :) = ad(j).firsths_dir45_eb;
    end
    gnames = cellfun(@num2str, num2cell(ad.hstx45),'UniformOutput',false);
    ylbl = 'probability head sweep is to left';
    for j = 1:nexp
        le{(2*j) - 1} = [po(j).legendEntry ' all hs'];
        le{(2*j)} = [po(j).legendEntry ' first hs'];
    end
    ttl = ('Head Swing Initiation Bias');
    bardat = hsbar';
    ebdat = hseb';
    hhh = barweb(bardat, ebdat, 0.8, gnames,ttl,  [],ylbl,[],[],le,[],'axis');
    for j = 1:nexp
        set(hhh.bars(2*j - 1), 'FaceColor', po(j).color, 'EdgeColor', po(j).color);
        set(hhh.bars(2*j), 'FaceColor', 'w', 'EdgeColor', po(j).color, 'LineWidth', 4);
    end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
    set(hhh.leg, 'FontSize', fontsize-2, 'FontName', font);
end


saveName =  'ReorientationRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    for j = 1:nexp
        errorbar (ad(j).txc, ad(j).reorate_thetabound, ad(j).reorate_thetabound_eb, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('reorientation rate (min$^{-1}$)');
    title ('Reorientation Rate vs. Direction of Forward Movement Relative To Boundary');
    if(showlegend), legend ({po.legendEntry}, 'Location', 'NorthEast','Interpreter','Latex'); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'ReorientationRateVsDistance';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    leg = {};
    for j = 1:nexp
        errorbar (10*ad(j).distx, ad(j).reorate_disttobound_todark, ad(j).reorate_disttobound_todark_eb, 'k--', 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
        errorbar (10*ad(j).distx, ad(j).reorate_disttobound_tolight, ad(j).reorate_disttobound_tolight_eb, 'k-', 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
        leg{2*j - 1} = [po(j).legendEntry ' towards dark'];
        leg{2*j} = [po(j).legendEntry ' towards light'];

    end
    xlim ([-5 5]);
    xlabel ('distance from boundary (mm)');
    ylabel ('reorientation rate (min$^{-1}$)');
    title ('Reorientation Rate vs. Distance of Head from Boundary');
    if(showlegend), legend ({po.legendEntry}, 'Location', 'NorthEast','Interpreter','Latex'); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end



saveName =  'ReorientationMagVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    for j = 1:nexp
        errorbar (ad(j).txc, ad(j).reosize_vs_ttb, ad(j).reosize_vs_ttb_eb, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('mean reorientation magnitude (deg)');
    title ('Size of Reorientation vs. Direction of Forward Movement Relative To Boundary');
    if(showlegend), legend ({po.legendEntry}, 'Location', 'NorthEast','Interpreter','Latex'); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'ReorientationDirVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    for j = 1:nexp
        errorbar (ad(j).txc, ad(j).reodir_vs_ttb, ad(j).reodir_vs_ttb_eb, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('mean reorientation (deg)');
    title ('Direction of Reorientation vs. Direction of Forward Movement Relative To Boundary');
    if(showlegend), legend ({po.legendEntry}, 'Location', 'NorthEast','Interpreter','Latex'); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end

saveName =  'PauseRateVsHeading';
if (isempty(whichGraphs) || any(strcmpi(saveName, whichGraphs)))
    
    nextfigure();
    hset(plotNumber).fignum = fignum;
    hset(plotNumber).saveName = saveName;

    for j = 1:nexp
        errorbar (ad(j).txc, ad(j).pauserate_thetabound, ad(j).pauserate_thetabound_eb, 'Color', po(j).color, 'Marker', po(j).marker,  'LineWidth', po(j).lineWidth, po(j).plotOptions{:}); hold on;
    end
    xlabel ('instantaneous run heading (degrees)');
    ylabel ('pausing rate (min$^{-1}$)');
    title ('Pausing Rate vs. Direction of Forward Movement Relative To Boundary');
    if(showlegend), legend ({po.legendEntry}, 'Location', 'NorthEast','Interpreter','Latex'); end
    emsmallen(gca, 'FontSize', fontsize, 'Font', font);
end


exts = {'.tiff', '.eps', '.ai', '.fig', '.jpg'}; ftype = {'tiff', 'eps2c', 'ai', 'fig', 'jpeg'};
if (~isempty(SaveDirectory))
    for j = 1:length(hset);
        figure(hset(j).fignum);
        for k = 1:length(exts)        
            fname = fullfile(SaveDirectory, [hset(j).saveName exts{k}]);
            saveas(hset(j).fignum, fname, ftype{k});
        end
    end
end
function nextfigure()
fignum = evalin('caller', 'fignum+1;'); assignin('caller', 'fignum', fignum);
forprinting = evalin('caller', 'forprinting;');
if(forprinting)
    figureForPrinting(fignum);
else
    figure(fignum); clf(fignum);
end
set(fignum, 'Color', 'w');  