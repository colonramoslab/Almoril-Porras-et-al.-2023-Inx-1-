function eset = fromDirectories(varargin)
%function eset = fromDirecotries()
%function eset = fromDirectories(dirname)
%function eset = fromDirectories(dirlist)
%
%loads set of experiments from disk
%usage: eset = MartaExperimentSet.fromDirectories()
%              prompts user to select a directory, then loads all .spine
%              files from subdirectories
%       eset = MartaExperimentSet.fromDirectories(dirname)
%              loads all subdirectories of dirname that have .spine files
%       eset = MartaExperimentSet.fromDirectories({dir1, dir2, dir3})
%              loads experiments from dir1, dir2, dir3
if (isempty(varargin))
    dirlist = MartaExperimentSet.findSpineDirectories();
else
    if (~iscell(varargin{1}))
        dirlist = MartaExperimentSet.findSpineDirectories('basedir', varargin{1});
    else
        dirlist = varargin{1};
    end
end

for j = 1:length(dirlist)
    expt(j) = MartaExperiment.fromFile(dirlist{j});
end

eset = MartaExperimentSet;
eset.expt = expt;