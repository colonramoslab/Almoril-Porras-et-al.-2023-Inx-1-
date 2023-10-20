function [fg,bak] = makeForeground (fstub, inds, extension)
%function [fg,bak] = makeForeground (fstub, inds)

existsAndDefault('extension', '.jpg');
fg = imread([fstub num2str(inds(1)) extension]);
bak = fg;
for j = 2:length(inds)
    im = imread([fstub num2str(inds(j)) extension]);
    bak = min(bak, im);
    fg = max(fg, im);
end

fg = fg-bak;