function drawAllowedTransitions(sm)
%function drawAllowedTransitions(sm)

cla
draw_graph(sm.allowedTransitions, {sm.segmentationClusters.name});
hold off;