function gleanInfoFromDirectories(basedir)

p = genpath(basedir);
tok = regexp(p, '([^;]*);', 'tokens');
for j = 1:length(tok)
    dirs{j} = tok{j}{1};
end
whos dirs
for j = 1:length(dirs)
    disp (['processing ' dirs{j}]); ts1 = tic;
    d = dir(fullfile(dirs{j}, 'directoryInfo.mat'));
    %if (isempty(d))

        disp (['processing ' dirs{j}]); ts1 = tic;
        directoryInfo = getDirectoryInformation(dirs{j});
        if (~isempty(directoryInfo))
            save(fullfile(dirs{j},  'directoryInfo.mat'), 'directoryInfo');
        end
        toc (ts1);


        clear directoryInfo
   % end
end

