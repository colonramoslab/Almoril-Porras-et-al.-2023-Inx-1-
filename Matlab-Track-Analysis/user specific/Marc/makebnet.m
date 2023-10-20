names = {'run', 'kink', 'curl up', 'omega turn', 'un curl', 'reverse',  'going backwards', 'second reversal', 'pause'};
ln = length(names);
allowed = zeros(ln);
%run -> run, kink, curl up, reverse, pause
allowed(1,[1 2 3 6 9]) = 1;

%kink -> run, kink
allowed(2,[1 2]) = 1;

%curl up -> curl up, omega turn
allowed(3, [3 4]) = 1;

%omega turn -> ot, un curl
allowed(4, [4 5]) = 1;

%un curl - > uc, run
allowed(5, [1 5]) = 1;

%reverse -> reverse, going backwards
allowed(6, [6 7]) = 1;

%re-reverse -> re-reverse, run
allowed(8, [1 8]) = 1;

%going backwards -> going backwards, re-reverse, curl up
allowed(7, [3 7 8]) = 1;
%pause -> pause, run
allowed(9, [1 9]) = 1;

figure(1); clf(1);
draw_graph(allowed, names);

