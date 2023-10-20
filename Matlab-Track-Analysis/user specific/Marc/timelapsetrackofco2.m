%mytrack.expt.fn = 'G:\maggot data\Extracted Files\CO2 50 mL in 2 L
%air\20101021_CS1_tracks.bin'
%
%mytrack.locInFile = 65039784
cla
if (~exist('im','var'))
    [im,x,y] = trackTimeLapse(mytrack); 
end
pcolor(x,y,double(im)); shading flat; colormap gray; axis equal; axis tight
rh = rectangle ('Position', [0 700 157 15.7], 'FaceColor', 'w', 'LineStyle', 'none');
th = text(157/2, 716, '1 cm', 'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', 'Arial', 'FontWeight', 'Bold', 'FontSize', 24);
set(gca, 'XTick', [], 'YTick', []);
xlc = mean(get (gca, 'XLim'));
[nx,ny] = normalizedCoordinates(xlc + [50 -50], 700 + [0 0]);
%tah=annotation('textarrow',nx,ny,'String','Direction of Travel');%, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', 'Arial', 'FontWeight', 'Bold');