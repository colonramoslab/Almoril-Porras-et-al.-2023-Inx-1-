function fastMaggotDraw (pt,camcalinfo)
%function fastMaggotDraw (pt)
%hopefully will speed up drawing a maggot track point by using old-style
% (not oop) language

if (isempty(pt.imData))
    return;
end

x = double(pt.imOffset(1)-1) + (1:size(pt.imData,2));
y = double(pt.imOffset(2)-1) + (1:size(pt.imData,1));
cm = gray(64);
pcolor (x,y,double(pt.imData)); shading flat; colormap(cm); 
axis equal
axis tight

if (~exist('camcalinfo', 'var'))
    camcalinfo = [];
end


h = realPtsToCamera(pt.head, camcalinfo);
m = realPtsToCamera(pt.mid, camcalinfo);
t = realPtsToCamera(pt.tail, camcalinfo);
c = realPtsToCamera(pt.contour, camcalinfo);
c(:,end+1) = c(:,1); %complete contour

hold on;
plot (h(1),h(2),'g*', t(1),t(2),'rh', [t(1),m(1),h(1)], [t(2),m(2),h(2)],'y-', c(1,:), c(2,:), 'r-');
hold off;
