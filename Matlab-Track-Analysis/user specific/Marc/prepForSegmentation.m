gq = GlobalQuantity();
gq.fieldname = 'lrdtheta';
gq.xField = 'theta';
gq.xData = 5; %smoothing time in seconds
gq.derivationMethod = @(xin, xData, yData) deriv(unwrap(xin), xData(1));

[cryo.expt.globalQuantity] = deal([]);
cryo.executeExperimentFunction('addGlobalQuantity', gq);
%%
tc = cryo.gatherSubField('sharpTurn', 'typeCode','expandToInds',true);
dtst = cryo.gatherFromSubField('sharpTurn','deltatheta');
ddst = cryo.gatherFromSubField('sharpTurn', 'ddtheta');
spst = cryo.gatherFromSubField('sharpTurn', 'speed');
scst = cryo.gatherFromSubField('sharpTurn', 'covRatio');
lrdtst = cryo.gatherFromSubField('sharpTurn', 'lrdtheta');


dtr = cryo.gatherFromSubField('run','deltatheta');
ddr = cryo.gatherFromSubField('run','ddtheta');
spr = cryo.gatherFromSubField('run', 'speed');
scr = cryo.gatherFromSubField('run', 'covRatio');
lrdtr = cryo.gatherFromSubField('run', 'lrdtheta');

clusterNames = {'omega turns', 'reversals', 'blips', 'sudden turns', 'runs'};
typeCode = [-1, 1, 0, -2, 0];
datafields = {'covRatio', 'deltatheta', 'speed', 'lrdtheta', 'ddtheta'};
operations = {@(x) x.^2, @(x) abs(x), @(x) x, @(x) abs(x), @(x) abs(x)};
clear sc;
for j = 1:length(clusterNames)
    sc(j) = SegmentationCluster(); %#ok<SAGROW>
    sc(j).name = clusterNames{j};
    sc(j).datafields = datafields;
    sc(j).operation = operations;
    sc(j).isSharpTurn = true;
    sc(j).typeCode = typeCode(j);
    sc(j).priorProbability = 1;
end
%%
data = [scst.^2; abs(dtst); abs(spst); abs(lrdtst); abs(ddst)];
for j = 1:3
    sc(j).clustMean = mean(data(:,tc == sc(j).typeCode),2);
    sc(j).clustCov = cov(data(:,tc == sc(j).typeCode)');
end
%%
data = [scr.^2; abs(dtr); abs(spr); abs(lrdtr); abs(ddr)];
%inds = all(isfinite(data), 1);

sc(4).clustMean = mean(data(:,all(isfinite(data), 1) & abs(lrdtr) > deg2rad(10)),2);
sc(4).clustCov = cov(data(:,all(isfinite(data), 1) & abs(lrdtr) > deg2rad(10))');

sc(5).clustMean = mean(data(:,all(isfinite(data), 1) & abs(lrdtr) < deg2rad(10)),2);
sc(5).clustCov = cov(data(:,all(isfinite(data), 1) & abs(lrdtr) < deg2rad(10))');

%%
figure(1); clf
track = cryo.expt(1).track(1);
[~,mostlikely] = max(sc.posteriorProbabilities(track));
c = 'rgbcmyk';
track.plotPath('sloc', 'b-', 'inds', mostlikely == 5); hold on
for j = 1:4
    track.plotPath('sloc', [c(j) '.'], 'inds', mostlikely == j);
    hold on
end
legend({sc([5 1:4]).name})
    
