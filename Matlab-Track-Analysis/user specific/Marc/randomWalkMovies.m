res = [768 1024];
avifn = {'C:\Documents and Settings\Marc\My Documents\presentations\randomWalk.avi', ...
    'C:\Documents and Settings\Marc\My Documents\presentations\biasLength.avi', ...
     'C:\Documents and Settings\Marc\My Documents\presentations\biasMag.avi', ...
      'C:\Documents and Settings\Marc\My Documents\presentations\biasDir.avi'};%output movie name
movieCodec = 'Cinepak'; % choose from 'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE' or 'None'
figdir = 'C:\Documents and Settings\Marc\My Documents\presentations\';
if (~exist('s','var'))

    s(1) = BiasedRandomWalkSimulator(200,5000);
    s(2) = BiasedRandomWalkSimulator(200,5000);
    s(3) = BiasedRandomWalkSimulator(200,5000);
    s(4) = BiasedRandomWalkSimulator(200,5000);
    
    s(1).setParams(deg2rad(1), 0, 1/20, 0, pi/2, 0, 0); s(1).driftCoeff = 0.1;
    s(2).setParams(deg2rad(1), 0, 1/20, 0.2, pi/2, 0, 0); s(2).driftCoeff = 0.1;
    s(3).setParams(deg2rad(1), 0, 1/20, 0.2, pi/2, 0.2, 0); s(3).driftCoeff = 0.1;
    s(4).setParams(deg2rad(1), 0, 1/20, 0.2, pi/2, 0.2, 0.2); s(4).driftCoeff = 0.1;

    tic
    for j = 1:length(s)
        s(j).simulate;
        toc;
    end
end
clear l c sig;
for j = 1:length(s)
    l{j} = permute(s(j).loc,[1 3 2]);
    c{j} = squeeze(mean(l{j},1));
    sig{j} = squeeze(std(l{j},1));
end

existsAndDefault('makeMovies', false);
existsAndDefault('showMovies', true);
existsAndDefault('saveFigs', false);

figure(2); clf(2);
set(2,'Color', 'k','InvertHardCopy','off');
plot (rad2deg(s(2).thetaAxis), s(2).reoRate, 'w-', 'LineWidth', 4);
xlim([0 360]);
ylim([0 1.1*max(s(2).reoRate)]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w','YTick', [],'XTick', 0:60:360);
xlabel ('direction');
ylabel ('probability of reorientation');
embiggen();
if (saveFigs)
    saveas(gcf, [figdir 'probability of reorientation.tiff'], 'tiff');
end
plot (rad2deg(s(3).thetaAxis), rad2deg(s(3).reoMag), 'w-', 'LineWidth', 4);
xlim([0 360]);
ylim([0 135]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w','YTick', [0 45 90 135],'XTick', 0:60:360);
xlabel ('direction');
ylabel ('magnitude of reorientation');
embiggen();
if saveFigs
    saveas(gcf, [figdir 'magnitude of reorientation.tiff'], 'tiff');
end
col = 'crgwymb';
tx = 1:length(c{1});
plot (tx, c{1}(:,1), col(1), tx, c{2}(:,1), col(2), tx, c{3}(:,1), col(3), tx, c{4}(:,1), col(4), ...
    'LineWidth', 4);
%{
hold on
for j = 1:4
    plot (tx, c{j}(:,1) + sig{j}(:,1), [col(j) ':'],  tx, c{j}(:,1) - sig{j}(:,1), [col(j) '--'], 'LineWidth', 2);
end
plot (tx, c{1}(:,1), col(1), tx, c{2}(:,1), col(2), tx, c{3}(:,1), col(3), tx, c{4}(:,1), col(4), ...
    'LineWidth', 4);
hold off
    %}
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w','YTick', [],'XTickLabel',[]);
leg = {'Random Walk', 'Biased Run Length', 'Biased Run Length + Reorientation Magnitude', 'Biased Run Length, Reorientation Mag, & Reorientation Direction'};
set(legend('Random Walk', 'Run Len', 'Run Len, Reo Mag', 'Run Len, Reo Mag + Dir', 'Location', 'NorthWest'),'FontSize',12,'TextColor', 'w');
        
xlabel ('time'); ylabel ('center of mass position');
if saveFigs
    saveas(gcf, [figdir 'center of mass travel.tiff'], 'tiff');
end



if (~makeMovies && ~showMovies)
    return;
end


figure(1); clf(1)
imshow(zeros(res));
ax = gca;
u = get(ax, 'Units');
set(ax, 'Units', 'Pixels');
avirect = get(ax, 'Position');
set(ax, 'Units', u);

axis(ax, [-res(1) res(1) -res(2) res(2)]); axis(ax, 'equal'), hold(ax,'off');
xl = get(ax, 'XLim');
yl = get(ax, 'YLim');
%col = 'brgmcyk';


for n = 1:4
    if makeMovies
        aviobj = avifile(avifn{n}, 'Quality', 90, 'Compression', movieCodec);
    else
        aviobj = [];
    end
    for j = 1:10:3000
        hold (ax, 'off');
        clear h;
        for k = 1:n
            if (k == n)
                ms = 16;
            else
                ms = 8;
            end
            plot (ax, l{k}(:,j,1), l{k}(:,j,2), [col(k) '.'], c{k}(1:j,1), c{k}(1:j,2), [col(k) '-'], 'LineWidth', 4, 'MarkerSize', ms);
            set(ax, 'XTick', [], 'YTick', [],'Color', 'k');
            hold (ax, 'on');
            h(k) = ellipse(3*sig{k}(j,1), 3*sig{k}(j,2), 0, c{k}(j,1), c{k}(j,2), col(k));
            set (h(k),'LineWidth', 3);
        end
        hold off
        set(legend(h, leg{1:n}, 'Location', 'NorthEast'),'FontSize',16,'TextColor', 'w');
        axis(ax, [xl yl]);
        set(ax, 'XTick', [], 'YTick', []);
        if makeMovies
            F = getframe(gcf, avirect);
            aviobj = addframe(aviobj, F);
        else
            pause(0.05);
        end
    end
    if (makeMovies)
        aviobj = close(aviobj);
    end
end
makeMovies = false;