function setupDirectories (username)
%function setupDirectories (username)
%adds paths to the matlab path so everything will run correctly
%(@Class directories do not need to be added to the path)
%if a username is given, adds user specific\username and subdirectories
%to the path
addpath (pwd);
addpath (fullfile(pwd, 'basic routines'));
addpath (fullfile(pwd,'example scripts'));
addpath (fullfile(pwd,'useful extra classes'));
addpath (fullfile(pwd,'MartaAnalysis'));
addpath (genpath(fullfile(pwd, 'AndyAnalysis')));
addpath (genpath(fullfile(pwd,'utility functions')));
addpath (genpath(fullfile(pwd,'guis')));
addpath (fullfile(pwd,'Simulation Classes'));
if (exist ('username', 'var') && ~isempty(username))
    addpath (genpath(fullfile(pwd, 'user specific', username, '')));
end
s = javaclasspath('-static');
for j = 1:length(s)
    [~,n{j}] = fileparts(s{j});
end

snyml = fullfile(pwd,'snakeyaml-1.6.jar');
if ~any(strcmpi('snakeyaml-1.6', n))
    javaaddpath (snyml);
end

%addpath (fullfile(pwd, 'bnt'));
if (exist ('genpathKPM', 'file'))
    addpath (genpathKPM (fullfile(pwd, 'bnt')));
end