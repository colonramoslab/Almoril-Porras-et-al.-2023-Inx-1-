
font = 'Arial';
fontsize = 10;

set(0,'DefaultAxesFontSize', fontsize);
set(0,'DefaultAxesFontName', font);

set(0, 'DefaultTextInterpreter', 'Latex');
figure(1); clf(1);
co2 = [4.21 4.02 3.93 3.82; 3.04 3.18 3.23 3.21; 2.51 2.49 2.53 2.51; 1.88 1.89 1.88 1.88; 1.3 1.25 1.23 1.26; 0.71 0.67 0.69 0.76; 0.26 0.42 0.51 0.56]';


x = 35*(3:-1:-3);
y = (0.5 + (0:3)*2.5)*25.4;

p = polyfit(x(2:end-1), mean(co2(:,2:end-1)),1);
linfit = p(1).*x + p(2);

h = plot (x/10,co2','LineWidth',3,'Marker','s','MarkerSize',10); hold on;
plot (x/10, linfit, 'k--', 'LineWidth', 3);
markers = 's^do';
for j = 1:4
    set(h(j),'Marker',markers(j));
end
for j = 1:length(y)
    leg{j} = [num2str(y(j)/10,'%.1f') ' cm from inlet'];
end
legend(leg,'Location', 'best')
xlabel ('distance from center (cm)');
ylabel ('carbon dioxide concentration ($\%$)');
title ('CO2 concentration vs. position in gradient apparatus')
%embiggen(gca,14);

figure(2); clf(2);
surf(x/10,y/10,co2); shading interp; hold on
[xx,yy] = meshgrid(x/10,y/10);
plot3(xx(:), yy(:), co2(:), 'r.', 'MarkerSize', 20);
axis ij;
xlabel ('distance from center (cm)');
ylabel ('distance from inlet (cm)');
zlabel ('co2 concentration ($\%$)');
hold off

figure(3); clf(3);
devco2 = co2-p(1).*repmat(x,[length(y) 1]) - p(2);
h = plot (x/10,devco2','LineWidth',3,'Marker','s','MarkerSize',10);
markers = 's^do';
for j = 1:4
    set(h(j),'Marker',markers(j));
end
for j = 1:length(y)
    leg{j} = [num2str(y(j)/10,'%.1f') ' cm from inlet'];
end
xlabel ('distance from center (cm)');
ylabel ('$\Delta$ carbon dioxide concentration ($\%$)');
title ('Deviation from linearity vs. position in gradient apparatus')
ylim([-0.6 0.6]);
legend(leg,'Location', 'best')

figure(4); clf(4);
w = 16.4;
h = 12.3;

pcolor(x/10,y/10, devco2/(2.5) * 100); shading interp; hold on
my = mean(y)/10;
hhh = plot ([-w/2 w/2 w/2 -w/2 -w/2], my + [-h/2 -h/2 h/2 h/2 -h/2],'r--','LineWidth', 3);
legend (hhh, 'Observation Region');
hold off
colorbar vert
axis equal
axis ij;
caxis ([-5 25]);
xlabel ('distance from center (cm)');
ylabel ('distance from inlet (cm)');
%clabel ('deviation from linearity ($\%$ of mean concentration)');

existsAndDefault ('savefiles', false);
if (savefiles)
    exts = {'.tiff', '.eps', '.ai', '.fig', '.jpg', '.pdf'}; ftype = {'tiff', 'eps2c', 'ai', 'fig', 'jpeg', 'pdf'};
    sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 calibrations';
    sn = {'co2 vs pos', 'co2 vs pos 2D', 'co2 dev from lin', 'co2 dev from lin 2D'};
    s = warning('off');    
    for j = 1:4
        figure(j);
        set(j, 'InvertHardcopy', 'off', 'Color', 'w');
       for k = 1:length(exts)        
            fname = fullfile(sd, [sn{j} exts{k}]);
            saveas(j, fname, ftype{k});
        end
    end
end