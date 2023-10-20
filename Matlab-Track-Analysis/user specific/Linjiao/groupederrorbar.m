function groupederrorbar(b, errdata)

h = bar('v6',b,'grouped');
%h = bar(b,'grouped'); If you are using MATLAB 6.5 (R13)
xdata = get(h,'XData');
sizz = size(b);

%determine the number of bars and groups
NumGroups = sizz(1);
SizeGroups = sizz(2);
NumBars = SizeGroups * NumGroups;

% Use the Indices of Non Zero Y values to get both X values
% for each bar. xb becomes a 2 by NumBars matrix of the X values.
INZY = [1 3];
xb = [];

for i = 1:SizeGroups
for j = 1:NumGroups
xb = [xb xdata{i}(INZY, j)];
end
end

%find the center X value of each bar.
for i = 1:NumBars
centerX(i) = (xb(1,i) + xb(2,i))/2;
end

% To place the error bars - use the following:
hold on;
%eh = errorbar(centerX,b,errdata); If you are using MATLAB 6.5 (R13)
eh = errorbar('v6',centerX,b,errdata);

set(eh(1),'linewidth',2); % This changes the thickness of the errorbars
set(eh(1),'color','g'); % This changes the color of the errorbars
set(eh(2),'linestyle','none'); % This removes the connecting line