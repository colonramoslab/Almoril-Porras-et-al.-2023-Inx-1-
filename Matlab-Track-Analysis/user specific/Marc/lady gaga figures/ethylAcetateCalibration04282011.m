xc_inches = 1.8662 + 1.378*(0:6) - 6; %x centers in inches
xc_cm = sort ([xc_inches - 0.2, xc_inches + 0.2] * 2.54);
yc_inches = 2.15 + 2.5*(0:3);
yc_cm = yc_inches * 2.54;
%readings in mv, from low to high

row{1} = [83.7  96.3  143.5 164.1 227 247 291 314 370 388 423 430 540 567];
row{2} = [94.1  106.3 154.5 176.8 236 260 307 330 383 399 435 456 525 555];
row{3} = [98.7  108.9 151.5 173.6 229 251 299 317 375 395 449 459 531 552];
row{4} = [101.2 111.6 151   172.3 224 241 294 306 360 380 431 NaN 504 519];

%row{4}(~isfinite(row(4))) = 
offset = 49.87;
mv_ppm = 0.96;
response_factor = 4.2;

outread = [300 325 310 305];
cf = mean(outread - offset)./(outread - offset);


readings = (cell2mat(row') - offset)/mv_ppm*response_factor;
cf = repmat(cf', 1, size(readings, 2));
figure(1); clf(1);

xc = xc_cm;
yc = yc_cm;

[xx,yy] = meshgrid(xc, yc);
F = TriScatteredInterp(xx(isfinite(readings(:))), yy(isfinite(readings(:))), readings(isfinite(readings(:))));
ireadings = reshape(F(xx(:), yy(:)), size(readings));


p = polyfit(xc(3:end-3), mean(readings(:,3:end-3)),1);
linfit = p(1).*xc + p(2);

h = plot (xc,readings,'LineWidth',3,'Marker','s','MarkerSize',10); hold on;
plot (xc, linfit, 'k--', 'LineWidth', 3);
markers = 's^do';
for j = 1:4
    set(h(j),'Marker',markers(j));
end
for j = 1:length(yc)
    leg{j} = [num2str(yc(j),'%.1f') ' cm from inlet'];
end
legend(leg{:}, 'linear fit - 90 ppm/cm','Location', 'best')
xlabel ('distance from center (cm)');
ylabel ('ethyl acetate concentration (ppm)');
title ('EtAc concentration vs. position in gradient apparatus')
%embiggen(gca,14);

figure(2); clf(2);
surf(xc,yc,ireadings); shading interp; hold on
[xx,yy] = meshgrid(xc,yc);
plot3(xx(:), yy(:), readings(:), 'r.', 'MarkerSize', 20);
axis ij;
xlabel ('distance from center (cm)');
ylabel ('distance from inlet (cm)');
zlabel ('co2 concentration ($\%$)');
hold off

figure(3); clf(3);
devEtAc = ireadings-p(1).*repmat(xc,[length(yc) 1]) - p(2);
h = plot (xc,devEtAc,'LineWidth',3,'Marker','s','MarkerSize',10);
markers = 's^do';
for j = 1:4
    set(h(j),'Marker',markers(j));
end
xlabel ('distance from center (cm)');
ylabel ('$\Delta$ EtAc concentration (ppm)');
title ('Deviation from linearity vs. position in gradient apparatus')
%ylim([-0.6 0.6]);
legend(leg,'Location', 'best')

figure(4); clf(4);
w = 16.4;
h = 12.3;

pcolor(xc,yc, devEtAc/mean(readings(isfinite(readings)))*100); shading interp; hold on
my = mean(yc);
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
    sd = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate calibrations';
    sn = {'EtAc vs pos', 'EtAc vs pos 2D', 'EtAc dev from lin', 'EtAc dev from lin 2D'};
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