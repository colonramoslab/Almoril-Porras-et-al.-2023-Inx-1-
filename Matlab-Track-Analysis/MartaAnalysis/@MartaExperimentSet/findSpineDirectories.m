function dirlist = findSpineDirectories(varargin)
%function dirlist = findSpineDirectories(varargin)

basedir = [];
dirlist = {};
verbose = false;
varargin = assignApplicable(varargin);

if (isempty(basedir))
    basedir = uigetdir([], 'Select a directory to search for spine subdirectories');
end
if (verbose)
    disp (['starting to parse ' basedir ', dirlist has ' num2str(length(dirlist)) ' entries']);
end
d = dir(fullfile(basedir, '*.spine'));
if ~isempty(d)
    if (verbose)
        disp (['dirlist has ' num2str(length(dirlist)) ' entries']);
        disp (['adding ' basedir ' to list']);
    end
    dirlist = [dirlist, basedir];
    if (verbose)
        disp (['dirlist has ' num2str(length(dirlist)) ' entries']);
    end
    return;
end

if (verbose)
    disp(['scanning ' basedir ' for subdirectories']);
end
%d = dir(fullfile(basedir,'*.'));
d = dir(basedir);
if (~isempty(d))
    if (strcmpi(d(1).name, '.') || strcmpi(d(1).name, '..'))
        d = d(2:end);
    end
end
if (~isempty(d))
    if (strcmpi(d(1).name, '.') || strcmpi(d(1).name, '..'))
        d = d(2:end);
    end
end
if (~isempty(d))
    d = d([d.isdir]);
end
if (isempty(d))
    return;
end
if (verbose)
    disp (['found ' num2str(length(d)) ' sub directories']);
end
for j = 1:length(d)
    if (~strcmp(d(j).name, '.') && ~strcmp(d(j).name, '..'))
        temp{j} = MartaExperimentSet.findSpineDirectories('basedir', fullfile(basedir, d(j).name), 'dirlist', dirlist);      
    end
end
dirlist = [dirlist, temp{:}];
    