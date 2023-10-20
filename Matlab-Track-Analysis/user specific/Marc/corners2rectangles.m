function [rectset1,rectset2] = corners2rectangles(xc,yc,ncols)
%function [rectset1,rectset2] = corners2rectangles(xc,yc,ncols)
%
%turns checkerboard corners into sets of rectangles;  up to you to figure
%out which is light and which is dark
%each set of rectangles follows the rules of structvar
%
%about structvar
%the first rule of structvar is you do not talk about structvar
%the second rule of structvar is you do not talk about structvar
%the third rule of structvar is that it is a Nx2 array of points,
%specifying rectangles
%the fourth rule of structvar is the lower left (in image coordinates; upper left in xy coords)
%corner comes first, then the other 3 points are specified in
%counterclockwise (in image coords; clockwise in xy coords) order
%the fifth rule of structvar is you do not talk about structvar
%

%sort points by row from high y to low y, then within rows from low x to high x
[yc,I] = sort(yc,'descend');
xc = xc(I);

for j = 1:(length(xc)/ncols)
    inds = ncols*(j-1) + (1:ncols);
    x = xc(inds);
    y = yc(inds);
    [x,I] = sort(x,'ascend');
    row(j).x = x;
    row(j).y = y(I);    
end
xc = [row.x];
yc = [row.y];

%determine the indices of the lower left corners of the two sets
llind1 = [];
llind2 = [];
for j = 1:(length(row) - 1)
    stagger = mod(j,2);
    llind1 = [llind1 ncols*j  + ((1+ stagger):2:(ncols-1))];
    llind2 = [llind2 ncols*j  + ((2-stagger):2:(ncols-1))];
end

%expand to take the rest of the corners in correct order
allinds1 = zeros([1 length(llind1) * 4]);
allinds1(1:4:end) = llind1;
allinds1(2:4:end) = llind1 + 1;
allinds1(3:4:end) = llind1 + 1 - ncols;
allinds1(4:4:end) = llind1 - ncols;

allinds2 = zeros([1 length(llind2) * 4]);
allinds2(1:4:end) = llind2;
allinds2(2:4:end) = llind2 + 1;
allinds2(3:4:end) = llind2 + 1 - ncols;
allinds2(4:4:end) = llind2 - ncols;

%assign to rectangles, keeping dimension same as structvar
if (size(xc,2) > 1)
    xc = xc';
end
if (size(yc,2) > 1)
    yc = yc';
end
rectset1 = [xc(allinds1),yc(allinds1)];
rectset2 = [xc(allinds2),yc(allinds2)];
    