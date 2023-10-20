basedir = 'E:\worm thermotaxis bin files';

p = genpath(basedir);
tok = regexp(p, '([^;]*);', 'tokens');
for j = 1:length(tok)
    dirs{j} = tok{j}{1};
end
if (~exist('di', 'var'))
    for j = 1:length(dirs)
        fname = fullfile(dirs{j}, 'directoryInfo.mat');
        d = dir(fname);
        if (~isempty(d))
            load (fname);
            di{j} = directoryInfo;
            
        end
    end
end

ur = [];
ll = [];
sn = [];
mt = [];
dicat = [];
for j = 1:length(di)
    if (~isempty(di{j}))
        dicat = [dicat; di{j}];
        sn = [sn [di{j}.setupNumber]];
        mt = [mt [di{j}.metrics]];
    end
end

ll = [mt.lowerleft];
ur = [mt.upperright];
com = [mt.centerOfMass];
x0 = [dicat.x0];
y0 = [dicat.y0];
plot (ll(1, sn == 1), ll(2, sn == 1), 'r.', ll(1, sn == 2), ll(2, sn == 2), 'g.', x0(sn == 0), y0(sn == 0), 'b.', ur(1, sn == 1), ur(2, sn == 1), 'r.', ur(1, sn == 2), ur(2, sn == 2), 'g.', x0(sn == 0), y0(sn == 0), 'b.');
hold on
plot (com(1, sn == 1), com(2, sn == 1), 'r.', com(1, sn == 2), com(2, sn == 2), 'g.', com(1, sn == 0), com(2, sn == 0), 'b.');
hold on
plot (ll(1, sn < 0), ll(2, sn < 0), 'k.', ur(1, sn < 0), ur(2, sn < 0), 'k.')
