function processAllLinjiaoToMatfiles(basedir)

p = genpath(basedir);
tok = regexp(p, '([^;]*);', 'tokens');
for j = 1:length(tok)
    dirs{j} = tok{j}{1};
end
for j = 1:length(dirs)
    d = dir(fullfile(dirs{j}, '*.bin'));
    if (isempty(d))
        continue;
    end
    
    
        disp (['processing ' dirs{j}]); ts1 = tic;
        processLinjiaoFilesToMatfiles(dirs{j});
        toc (ts1);
end

