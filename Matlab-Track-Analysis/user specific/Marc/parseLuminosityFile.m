function lum = parseLuminosityFile (fn)
%function lum = parseLuminosityFile (fn)

data = load(fn);
for j = 1:size(data,1)
    lum(j).box = data(j,1:4)';
    lum(j).cx = mean(lum(j).box([1 3]));
    lum(j).cy = mean(lum(j).box([2 4]));
    lum(j).proj = data(j,5:3:end)';
    lum(j).lux = data(j,6:3:end)';
    lum(j).cam = data(j,7:3:end)';
end

%nuke all 0s, something went wrong
valid = any([lum.lux] > 0);
lum = lum(valid);

for j = 1:length(lum)
    p = polyfit(lum(j).lux, lum(j).cam, 1);
    lum(j).m = p(1);
    lum(j).b = p(2);
end

for j = 1:length(lum)
    x0 = [0.01 2.35 20];
    lum(j).p2l = lsqcurvefit(@(x,xdata) x(1).*xdata.^(x(2)) + x(3), x0, lum(j).proj, lum(j).lux);
end