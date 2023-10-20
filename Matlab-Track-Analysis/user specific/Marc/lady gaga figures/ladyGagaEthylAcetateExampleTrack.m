basedirs = {'E:\larvalco2\Extracted\Ethyl Acetate 233ppm CS\50mL air in 1.5 L air\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 2.25ppm CS\40mL air in 2 L air\CS', ...
    'E:\larvalco2\Extracted\Ethyl Acetate masking experiment\25ppm masking 25ppm in vlaves\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 331ppm Gr63a\50mL air in 1 L air\Gr63a',...
    'E:\larvalco2\Extracted\Ethyl Acetate pure Spatial Gr63a\50mL co2 in 2 L air\Gr63a', 'E:\larvalco2\Extracted\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS',...
    'E:\larvalco2\Extracted\Ethyl Acetate Or83b1\23ppm at output spatial\Or83b1', 'E:\larvalco2\Extracted\ethyl acetate Or83b2 spatial\23ppm middle in 2L air\Or83b2'}; 

esetnames = {'cs_ea_0_460', 'cs_ea_0_5', 'cs_ea_25_50', 'gr63a_ea_0_660', 'gr63a_ea_50_2', 'cs_ea_0_48', 'or83b_ea_0_46','Or83b2_ea_0_46'};
ecl = ESetCleaner();
ecl.rpmCut = 1;
ecl.askFirst = false;
ecl.showFigsInReport = false;
        

if (~exist('eset', 'var'))
    j = 6;
    eset = ExperimentSet.fromMatFiles(fullfile (basedirs{j}, 'matfiles', esetnames{j}),1);
    ecl.getReport(eset);
    ecl.clean(eset);
end
if (~exist('t', 'var'))
    t = [eset.expt.track];
end
%tind = find([t.locInFile] == 412303112, 1, 'first');
% tind = find([t.locInFile] == 170864024, 1, 'first');
% 
% 
% if (isempty(tind))
%     disp ('track not found');
%     return;
% end
% 

tind = 10;
sh = 'interp';
ms = 8;
track = t(tind);

inds = [30   444];

track.setSegmentSpeeds;
track.segmentTrack;
[im, x, y] = trackTimeLapse(track, 'indsRange', inds);
thresh = 30;
imt = im > thresh;
imt = imdilate(imt, ones(7));
imt = imclose(imt, ones(7));
im = im.*imt;
%im(im < thresh) = 0;

figure(1); clf();
clear l tt;
%{
p = get(gcf, 'Position');

p(4) = p(3); %*9/8
if (p(2) + p(4) > 1024)
    p(2) = 1024 - p(4);
end

set(gcf, 'Position', p);
%}
for j = 1:3
    aax(j) = subplot(3,1,j);
end
p = get(aax(3), 'Position');
y0 = p(2) - 0.05;
p = get(aax(1), 'Position');
h = (p(2) + p(4) + 0.05 - y0)/3;
for j = 1:3
    p = get(aax(j), 'Position');
    p(2) = y0 + (3-j)*h;
    p(4) = 0.95*h;
    set (aax(j), 'Position', p);
end

axes(aax(1));

pcolor(x, y, double(im)); shading (sh); 

cm = gray(256);
cm = cm(end:-1:1,:);
colormap (cm);
hold on;
sl = track.expt.camcalinfo.camPtsFromRealPts(track.getDerivedQuantity('sloc'));
plot (sl(1,inds(1):inds(end)), sl(2,inds(1):inds(end)), 'b.', 'MarkerSize', ms/2);
plot (sl(1,inds(1):5:inds(end)), sl(2,inds(1):5:inds(end)), 'b.', 'MarkerSize', ms, 'LineWidth', 1.5); 

ratio = 4.86;
%xc = 1628;

%xlim(xc + [-1/2 1/2]*ratio*80);
%xlim([1420 1720]);
ylim([700 790]);
set(gca, 'DataAspectRatio', [1 1 1]);
axis(gca, 'fill');
%xc = 1674;
%yc = 724;
%r = 23;

xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
l(1) = plot ([xl(1) + 10, xl(1) + 10 + track.expt.camcalinfo.pixelsPerRealUnit/2], [yl(2) - 10, yl(2) - 10], 'k-', 'LineWidth', 5);
tt(1) = text ( xl(1) + 10 + track.expt.camcalinfo.pixelsPerRealUnit/4, yl(2) - 6, '{\bf 1 cm}','VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'latex');
xl = 11.5875 + [-0.15 0.15];
yl = 4.7727 + [-0.15 0.15];
pts = [xl(1) xl(1) xl(2) xl(2) xl(1); yl(1) yl(2) yl(2) yl(1) yl(1)];
pts = track.expt.camcalinfo.camPtsFromRealPts(pts);
pts(1, [1 2 5]) = min(pts(1, [1 2 5]));
pts(1, [3 4]) = max(pts(1, [3 4]));
pts(2, [1 4 5]) = min(pts(2, [1 4 5]));
pts(2, [2 3]) = max(pts(2, [2 3]));
%plot (xc + r*cos(linspace(0, 2*pi, 360)), yc + r*sin(linspace(0, 2*pi, 360)), 'r--', 'LineWidth', 2);
l(2) = plot (pts(1,:), pts(2,:), 'k-.', 'LineWidth', 1);

    
xl2 = 11.5408 + [-0.15 0.15];
yl2 = 5.0294 + [-0.15 0.15];
pts = [xl2(1) xl2(1) xl2(2) xl2(2) xl2(1); yl2(1) yl2(2) yl2(2) yl2(1) yl2(1)];
pts = track.expt.camcalinfo.camPtsFromRealPts(pts);
pts(1, [1 2 5]) = min(pts(1, [1 2 5]));
pts(1, [3 4]) = max(pts(1, [3 4]));
pts(2, [1 4 5]) = min(pts(2, [1 4 5]));
pts(2, [2 3]) = max(pts(2, [2 3]));
%plot (xc + r*cos(linspace(0, 2*pi, 360)), yc + r*sin(linspace(0, 2*pi, 360)), 'r--', 'LineWidth', 2);
l(3) = plot (pts(1,:), pts(2,:), 'r:', 'LineWidth', 1);


hold off
ax1 = gca;
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 1);

movinds = track.getDerivedQuantity('mapinterpedtopts', false,[296 304 312 320]);
axes(aax(2));
ax2 = gca();
pos = get(ax1, 'Position');
pbar = get(ax1, 'PlotBoxAspectRatio');
if (pos(3) > pos(4) * pbar(1)/pbar(2))
    w = pos(3) * pbar(1)/pbar(2);
else
    w = pos(3)/4;
end
x0 = pos(1) + pos(3)/2 - 2*w;
for j = 2:4
    ax2(j) = cloneaxes(ax2(1));
end

for j = 1:4
    p = get(ax2(1), 'Position');
    p(3) = w ;%* 0.98;
    p(1) = x0 + (j-1) * w;
    set(ax2(j), 'Position', p);
end
for j = 1:4
    track.pt(movinds(j)).drawTrackImage(track.expt.camcalinfo, 'fid', track.expt.fid, 'Axes', ax2(j), 'pretty', true, 'drawHeadArrow', false, 'spineColor', 'y.-', 'contourWidth', 1); 
    hold(ax2(j), 'on');
    track.plotPath('sloc', 'b.', 'MarkerSize', ms/2, 'Axes', ax2(j), 'inds', inds(1):inds(end));
    track.plotPath('sloc', 'b.', 'MarkerSize', ms, 'Axes', ax2(j), 'inds', inds(1):5:inds(end));
    shading(ax2(j), sh);
    axis(ax2(j), 'fill');
    set(ax2(j), 'XTick', [], 'YTick', [], 'LineWidth', 1, 'XLim', xl2, 'YLim', yl2, 'DataAspectRatio', [1 1 1], 'Visible', 'off');
    %annotation('rectangle', plotboxpos(ax(j)),'LineWidth',1,'LineStyle','-', 'EdgeColor', 'k'); 
end
bb = plotboxpos(ax2(1));
pb = plotboxpos(ax2(4));

bb(3) = pb(1) + pb(3) - bb(1);
bb(4) = pb(2) + pb(4) - bb(2);
b(1) = annotation('rectangle',bb,'LineWidth',1,'LineStyle',':', 'EdgeColor', 'r'); 
for j = 1:3
    l = [l annotation('line', j*pb(3) + [bb(1) bb(1)], bb(2) + [0 bb(4)], 'LineWidth',1,'LineStyle',':', 'Color', 'r')]; 
end
et = [track.pt(movinds).et];
et0 = et(1);
et = et - et0;

for j = 1:4
    xltemp = get(ax2(j), 'XLim');
    yltemp = get(ax2(j), 'YLim');
    tt = [tt text(xltemp(1) + 0.01, yltemp(2) - 0.01, ['{\bf t = ' num2str(et(j), '%.1f') ' s}'],'VerticalAlignment', 'top', 'Parent', ax2(j), 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'latex')];
end

axes(aax(3))
movinds = track.getDerivedQuantity('mapinterpedtopts', false, [351, 362, 373, 384]);
ax = gca();
pos = get(ax1, 'Position');
pbar = get(ax1, 'PlotBoxAspectRatio');
if (pos(3) > pos(4) * pbar(1)/pbar(2))
    w = pos(3) * pbar(1)/pbar(2);
else
    w = pos(3)/4;
end
x0 = pos(1) + pos(3)/2 - 2*w;
for j = 2:4
    ax(j) = cloneaxes(ax(1));
end

for j = 1:4
    p = get(ax(1), 'Position');
    p(3) = w ;%* 0.98;
    p(1) = x0 + (j-1) * w;
    set(ax(j), 'Position', p);
end
for j = 1:4
    track.pt(movinds(j)).drawTrackImage(track.expt.camcalinfo, 'fid', track.expt.fid, 'Axes', ax(j), 'pretty', true, 'drawHeadArrow', false, 'spineColor', 'y.-', 'contourWidth', 1); 
    hold(ax(j), 'on');
    track.plotPath('sloc', 'b.', 'MarkerSize', ms/2, 'Axes', ax(j), 'inds', inds(1):inds(end));
    track.plotPath('sloc', 'b.', 'MarkerSize', ms, 'Axes', ax(j), 'inds', inds(1):5:inds(end));
    shading(ax(j), sh);
    axis(ax(j), 'fill');
    set(ax(j), 'XTick', [], 'YTick', [], 'LineWidth', 1, 'XLim', xl, 'YLim', yl, 'DataAspectRatio', [1 1 1], 'Visible', 'off');
    %annotation('rectangle', plotboxpos(ax(j)),'LineWidth',1,'LineStyle','-', 'EdgeColor', 'k'); 
end
bb = plotboxpos(ax(1));
pb = plotboxpos(ax(4));

bb(3) = pb(1) + pb(3) - bb(1);
bb(4) = pb(2) + pb(4) - bb(2);
b(2) = annotation('rectangle',bb,'LineWidth',1,'LineStyle','-.', 'EdgeColor', 'k'); 
for j = 1:3
    l = [l annotation('line', j*pb(3) + [bb(1) bb(1)], bb(2) + [0 bb(4)], 'LineWidth',1,'LineStyle','-.', 'Color', 'k')]; 
end
et = [track.pt(movinds).et];
et = et - et0;

for j = 1:4
    xl = get(ax(j), 'XLim');
    yl = get(ax(j), 'YLim');
    tt= [tt text(xl(1) + 0.01, yl(2) - 0.01, ['{\bf t = ' num2str(et(j), '%.1f') ' s}'],'VerticalAlignment', 'top', 'Parent', ax(j), 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'latex')];
end

colormap (cm);
set(gcf, 'Color', 'w');


exts = {'.tiff', '.eps', '.ai', '.fig', '.jpg', '.pdf'}; ftype = {'tiff', 'eps2c', 'ai', 'fig', 'jpeg', 'pdf'};
existsAndDefault('savefig', false);
SaveDirectory = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\cs ea 2 ppm per cm';
if (savefig)
    s = warning('off');   
    set(gcf, 'InvertHardcopy', 'off');
    for k = 1:length(exts)
        fname = fullfile(SaveDirectory, ['example track', exts{k}]);
        saveas(gcf, fname, ftype{k});
    end
    warning(s);

    colormap gray;
    set(gcf, 'Color', 'k');
    set([ax1, ax, ax2] , 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    set(l, 'Color', 'w');
    set(tt, 'Color', 'w');
    set(b, 'Color', 'w');
    s = warning('off');   
    
    set(gcf, 'InvertHardcopy', 'off');
    for k = 1:length(exts)
        fname = fullfile(SaveDirectory, ['example track inverted', exts{k}]);
        saveas(gcf, fname, ftype{k});
    end
    warning(s);

end



