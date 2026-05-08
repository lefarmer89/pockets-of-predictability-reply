function fp = saveFigureScript(funcName, ~, outDir)
%SAVEFIGURESCRIPT  Run a figure-rendering function and save the figure as
% a .eps file.
%
%   fp = saveFigureScript(funcName, [unused], outDir)
%
% funcName  e.g. 'figure1'. A function in figures/ that creates a MATLAB
%           figure (gcf valid after running).
% outDir    Output directory for the .eps file. Defaults to
%           'output/figures' relative to the package root.
%
% Returns the absolute path of the written .eps file.

if nargin < 3 || isempty(outDir); outDir = fullfile('output', 'figures'); end
if ~isfolder(outDir); mkdir(outDir); end

fp = fullfile(outDir, [funcName, '.eps']);

% Suppress the figure script's stdout (most don't write any, but a few
% have debug `disp` calls that shouldn't mix into displayResults output).
evalc(funcName);

if isempty(get(0, 'Children'))
    warning('saveFigureScript:NoFigure', ...
        '%s did not create a figure; skipping save.', funcName);
    fp = '';
    return
end

print(gcf, fp, '-depsc');
fprintf('=== %s  →  %s ===\n', funcName, fp);
end
