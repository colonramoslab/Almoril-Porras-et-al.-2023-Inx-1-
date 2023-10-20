function expt = fromFile (dirname)
%loads a set of tracks from a directory
%expt = MartaExperiment.fromFile (dirname);
%
%directory must contain .spine files
%may optionally contain .outline files

d = dir(fullfile(dirname, '*.spine'));
if (isempty(d))
    disp (['no .spine files found in ' dirname]);
    expt = [];
    return;
end
disp (['loading ' num2str(length(d)) ' tracks']);
ts = tic;
for j = 1:length(d)
    if (fix(4 * j / length(d)) > fix(4 * (j-1)/length(d)))
        disp ([num2str(toc(ts), 3), ' s elapsed, ' num2str(100 * j / length(d), 2) '% done']);
    end
    t(j) = MartaTrack.fromFile(fullfile(dirname, d(j).name));
end

expt = MartaExperiment;
expt.fname = dirname;
expt.dr = DerivationRules;
expt.dr.interpTime = 0.025;
expt.dr.smoothTime = 0.05;
expt.dr.derivTime = 0.025;

expt.so = MaggotSegmentOptions;
expt.so.speed_field = 'smoothSpeed';
expt.so.minRunTime = 1;
expt.so.minStopTime = 1;
expt.so.curv_cut = 8;

[t.dr] = deal([expt.dr]);
[t.so] = deal([expt.so]);
%[t.expt] = deal(expt);

expt.track = t;