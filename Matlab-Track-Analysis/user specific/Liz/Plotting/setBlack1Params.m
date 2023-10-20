% set line weights and font sizes
set(gca,'LineWidth',2);
set(gca,'FontSize',16);

% set colors
set(gcf,'Color',[0 0 0]);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);

% set overall figure size
set(gcf,'Position',[360 239 369+200 239+200]);

% set x and y axis label font size
xhand= get(gca, 'xlabel');
set(xhand, 'FontSize', 16);
yhand= get(gca, 'ylabel');
set(yhand, 'FontSize',16);
