npts = 1E5;

tic;
x = 1:npts;
toc

tic;
x2 = zeros([1 npts]);
for j = 1:npts
    x2(j) = j;
end
toc

tic;
for j = 1:npts
    x3(j) = j;
end
toc