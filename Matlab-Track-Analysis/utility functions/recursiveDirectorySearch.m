function fileNames = recursiveDirectorySearch(startingDirectory, fileSpecifier)
%function fileNames = recursiveDirectorySearch(startingDirectory, fileSpecifier)
%
%returns complete paths to all files in startingDirectory or subdirectories
%that satisfy fileSpecifier
%
%example
%startingDirectory = '\\labnas1\share\David\Extracted\Spatial'
%batfiles = recursiveDirectorySearch(startingDirectory, '*.bat')

p = genpath(startingDirectory);
inds = [0 strfind(p,';')];
for j = 1:(length(inds) - 1)
    directories{j} = p((inds(j)+1):(inds(j+1) - 1));
end

fnnum = 0;
for j = 1:length(directories)
    d = dir([directories{j} '\' fileSpecifier]);
    for k = 1:length(d)
        fnnum = fnnum+1;
        fileNames{fnnum} = [directories{j} '\' d(k).name];
    end
end
