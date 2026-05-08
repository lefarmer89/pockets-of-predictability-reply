function paths = setup_paths()
% SETUP_PATHS  Add package subfolders to the MATLAB path and return a
% struct of canonical directories. Idempotent. Returns paths regardless of
% whether they were already on the path.
%
%   paths = setup_paths()
%
% Fields of paths:
%   root           package root
%   subroutines    helper functions (organized into 7 topical subdirs:
%                  forecasting, diagnostics, portfolio, bootstrap,
%                  hyperparameters, io, utils)
%   tables         table*.m and tableA*.m
%   figures        figure*.m and figureA*.m
%   data           input data (predictors, calibrations)
%   csvSims        large simulation CSVs
%   results        cached pipeline outputs (.mat)
%   output         rendered tables/figures and logs
%   outputTables   final .tex
%   outputFigures  final .eps
%   outputLogs     diary logs

here = fileparts(mfilename('fullpath'));

addpath(here);
addpath(fullfile(here, 'subroutines'));
subDirs = {'forecasting','diagnostics','portfolio','bootstrap', ...
           'hyperparameters','io','utils'};
for k = 1:numel(subDirs)
    addpath(fullfile(here, 'subroutines', subDirs{k}));
end
addpath(fullfile(here, 'tables'));
addpath(fullfile(here, 'figures'));

paths = struct( ...
    'root',          here, ...
    'subroutines',   fullfile(here, 'subroutines'), ...
    'tables',        fullfile(here, 'tables'), ...
    'figures',       fullfile(here, 'figures'), ...
    'data',          fullfile(here, 'data'), ...
    'csvSims',       fullfile(here, 'data', 'csv_sims'), ...
    'results',       fullfile(here, 'results'), ...
    'output',        fullfile(here, 'output'), ...
    'outputTables',  fullfile(here, 'output', 'tables'), ...
    'outputFigures', fullfile(here, 'output', 'figures'), ...
    'outputLogs',    fullfile(here, 'output', 'logs') ...
);
end
