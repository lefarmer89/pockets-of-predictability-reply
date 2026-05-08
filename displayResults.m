function displayResults()
%DISPLAYRESULTS  Render every reply table and figure from cached .mat
% artifacts.
%
%   displayResults()
%
% For each table in `tables/`, this captures the LatexTableFull output
% and saves it to `output/tables/<name>.tex` (full table) plus
% `output/tables/<name>_data.tex` (data rows only, for `\input{}` from
% the reply manuscript). For each figure in `figures/`, the figure is
% rendered and saved to `output/figures/<name>.eps`.
%
% Entry-point routines (`dailyEmpirics`, `monthlyEmpirics`, etc.) write
% the cached `.mat` artifacts; this script only renders.

setup_paths();

outTables = fullfile('output', 'tables');
outFigs   = fullfile('output', 'figures');
if ~isfolder(outTables); mkdir(outTables); end
if ~isfolder(outFigs);   mkdir(outFigs);   end

% --- Tables ------------------------------------------------------------
tableScripts = {
    'table1',  'Table 1';
    'table2',  'Table 2';
    'table3',  'Table 3';
    'table4',  'Table 4';
    'table5',  'Table 5';
    'tableA1', 'Table A.1';
    'tableA2', 'Table A.2';
    'tableA3', 'Table A.3';
    'tableA4', 'Table A.4';
    'tableA5', 'Table A.5';
    'tableA6', 'Table A.6'};
for k = 1:size(tableScripts, 1)
    fprintf('\n%s\n', tableScripts{k, 2});
    saveTableScript(tableScripts{k, 1}, [], outTables);
end

% --- Multi-panel data splits ------------------------------------------
% Each multi-panel table emits a single _data.tex with N rows summed
% across panels. Split into per-panel files for `\input{}` between the
% panel's multicolumn header rows in the reply LaTeX.
%
% Each row: {scriptName, panelSizes, stripTrailingAmp}. Strip the
% trailing `&` only for tables whose tabular spec has no spacer column.
multiPanelSplits = {
    'table1',  [9 9],   false;
    'table2',  [2 3 3], true;
    'table3',  [9 10],  false;
    'table4',  [7 7 7], false;
    'table5',  [7 7 7], false;
    'tableA1', [9 9],   true;
    'tableA2', [9 9],   false;
    'tableA3', [9 9 9], false;
    'tableA4', [9 10],  false;
    'tableA5', [2 3 3], true;
    'tableA6', [9 10],  false;
};
for k = 1:size(multiPanelSplits, 1)
    splitDataPanels(multiPanelSplits{k, 1}, multiPanelSplits{k, 2}, ...
        multiPanelSplits{k, 3});
end

% --- Figures -----------------------------------------------------------
figureScripts = {'figure1', 'figure2', 'figureA1', 'figureA2', 'figureA3'};
for k = 1:numel(figureScripts)
    saveFigureScript(figureScripts{k}, [], outFigs);
end

fprintf('\n--- displayResults complete ---\n');
fprintf('  LaTeX tables: %s\n', outTables);
fprintf('  EPS figures:  %s\n', outFigs);
end
