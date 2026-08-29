%PROJECT_CONFIG_EXAMPLE Local configuration for photon-transport studies.
% Copy this file to project_config.m and edit the paths for your machine.
% project_config.m is ignored by Git to avoid publishing local paths.

PROJECT_ROOT = fileparts(fileparts(mfilename('fullpath')));
DATA_ROOT = fullfile(PROJECT_ROOT, 'data');
RESULTS_ROOT = fullfile(PROJECT_ROOT, 'results');

% Set this to the directory that contains the MCmatlab package.
MCMATLAB_ROOT = '';

if ~isempty(MCMATLAB_ROOT)
    addpath(MCMATLAB_ROOT);
end

if ~exist(DATA_ROOT, 'dir')
    mkdir(DATA_ROOT);
end

if ~exist(RESULTS_ROOT, 'dir')
    mkdir(RESULTS_ROOT);
end
