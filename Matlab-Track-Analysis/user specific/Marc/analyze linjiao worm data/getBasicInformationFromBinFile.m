function metrics = getBasicInformationFromBinFile (fname)

expt = Experiment.fromFile (fname, 'doNOTloadTIME', false, [], 500);

allloc = double(expt.gatherField('loc'));
valid = isfinite(allloc(1,:)) & isfinite(allloc(2,:));
allloc = allloc(:,valid);
metrics.lowerleft = min(allloc, [],2);
metrics.upperright = max(allloc, [], 2);

metrics.centerOfMass = mean(allloc, 2);

x = allloc(1,:);
y = allloc(2,:);

metrics.xpctls = percentile(x, 0.1:0.1:0.9);
metrics.ypctls = percentile(y, 0.1:0.1:0.9);

metrics.xaxis = min(x):10:max(x);
metrics.yaxis = min(y):10:max(y);
metrics.xhist = hist(x,metrics.xaxis);
metrics.yhist = hist(y,metrics.yaxis);

metrics.hist2d = makeIm(x,y, metrics.xaxis, metrics.yaxis);
