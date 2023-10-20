if (~exist('eset', 'var'))
    eset = ExperimentSet.fromMatFiles('E:\larvalco2\Extracted\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS\matfiles\cs_ea_2000x_dilution_50_2_spatial');
    ecl = ESetCleaner;
    ecl.rpmCut = 1;
    ecl.askFirst = false;
    ecl.clean(eset);
end

if (~exist('t', 'var'))
    t = [eset.expt.track];
end
font = 'Arial';
fontsize = 8;
figure(1); clf();
%t([t.startFrame] < 300 & [t.npts] > 1500).plotPath('displacement','k-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5); axis equal; hold on
t.plotPath('displacement','k-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5); axis equal; hold on

%cc = 'rgbcmy';
cc = {[0.8 0 0], [0 0.8 0], [0 0 0.8], [0 0.6 0.6], [0.6 0 0.6], [0.5 0.5 0.1], [0.3 0.1 0.6]}; 
cind = 0;
for j = [12 14 30 44 50 68 70]
    cind = mod(cind, length(cc)) + 1;
    t(j).plotPath('displacement', 'k-', 'Color', cc{cind}, 'LineWidth', 2.0);
end
xlim([-10 10]);

barloc = [9, 5.5];
yl = get(gca, 'YLim');
barloc(2) = yl(2) - 0.3;
plot (barloc(1) + [-0.5 0.5], barloc([2 2]), 'k-', 'LineWidth', 4);
text (barloc(1), barloc(2)-0.2, '1 cm', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontName', font, 'FontSize', fontsize);
plot (0, 0, 'r.', 'MarkerSize', 30);
%plot (0, 0, 'ro', 'MarkerSize', 20, 'LineWidth', 4);

hold off

set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2);

existsAndDefault('saveFig', false);
if (saveFig)
    SaveDirectory = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\cs ea 2 ppm per cm';
    exts = {'.tiff', '.eps', '.ai', '.fig', '.jpg', '.pdf'}; ftype = {'tiff', 'epsc2', 'ai', 'fig', 'jpeg', 'pdf'};
    s = warning('off');

    figure(1);
    set(1, 'InvertHardcopy', 'off', 'Color', 'w');
    for k = 1:length(exts)
        fname = fullfile(SaveDirectory, ['ea 2ppm PER cm tracks' exts{k}]);
        saveas(1, fname, ftype{k});
    end
 
    warning(s);
end

%{
figure(2); clf();

tind = 26;
t(tind).setSegmentSpeeds;
t(tind).segmentTrack;
[im, x, y] = trackTimeLapse(t(tind));
%[xx,yy] = meshgrid(x,y);
%cpts = t(tind).expt.camcalinfo.realPtsFromCamPts([xx(:) yy(:)]');
%pcolor (reshape(cpts(1,:),size(xx)),reshape(cpts(2,:),size(yy)), im); shading flat; colormap gray
pcolor (x, y, im); shading interp; colormap gray
sl = t(tind).expt.camcalinfo.camPtsFromRealPts(t(tind).getDerivedQuantity('sloc'));
hold on;
plot (sl(1,:), sl(2,:), 'r-'); 
hold off
%t(tind).plotPath('iloc', 'r-', 'LineWidth', 1);
%}