function splitDataPanels(scriptName, panelSizes, stripTrailingAmp)
% SPLITDATAPANELS  Split a rendered _data.tex into per-panel files for
% multi-panel tables, so the reply LaTeX can \input each panel between
% its own multicolumn header rows.
%
%   splitDataPanels(scriptName, panelSizes)
%   splitDataPanels(scriptName, panelSizes, stripTrailingAmp)
%
% scriptName        e.g. 'table1' (no path; assumed under output/tables/).
% panelSizes        vector of row counts per panel, summing to the total
%                   number of data rows in <scriptName>_data.tex.
% stripTrailingAmp  optional logical (default false). When true, strips
%                   the trailing empty `&` cell that LatexTableFull
%                   always emits. Use this for tables whose tabular
%                   spec exactly matches the visible field count
%                   (Tables 2, A.1, A.5). Default keeps the trailing `&`
%                   to match the spec's spacer column (Tables 1, 3, 4,
%                   A.2, A.3, A.4, A.6).
%
% Writes <scriptName>_panel1.tex, _panel2.tex, ... (1-indexed). Strips
% the trailing `\\` row terminator from each panel's last line so the
% reply LaTeX can place its own `\\` before \midrule / \bottomrule.

if nargin < 3 || isempty(stripTrailingAmp); stripTrailingAmp = false; end

outDir = fullfile('output', 'tables');
inFile = fullfile(outDir, [scriptName, '_data.tex']);
if ~isfile(inFile)
    error('splitDataPanels:noFile', 'data file not found: %s', inFile);
end

txt   = fileread(inFile);
lines = regexp(txt, '\r?\n', 'split');
% Drop blank trailing element from final newline.
while ~isempty(lines) && isempty(strtrim(lines{end}))
    lines(end) = [];
end
if numel(lines) ~= sum(panelSizes)
    error('splitDataPanels:rowMismatch', ...
        '%s has %d rows but panelSizes sums to %d', ...
        inFile, numel(lines), sum(panelSizes));
end

cursor = 0;
for p = 1:numel(panelSizes)
    n        = panelSizes(p);
    panelLines = lines(cursor + (1:n));
    cursor   = cursor + n;
    if stripTrailingAmp
        % Strip the trailing empty `&` cell from every row when the
        % tabular spec doesn't include a spacer column.
        for k = 1:numel(panelLines)
            panelLines{k} = regexprep(panelLines{k}, '\s*&\s*(\\\\)?\s*$', ' $1');
            panelLines{k} = regexprep(panelLines{k}, '\s+$', '');
        end
    end
    % Strip the trailing `\\` from this panel's last line — the reply
    % LaTeX provides the row terminator after \input{}.
    panelLines{end} = regexprep(panelLines{end}, '\s*\\\\\s*$', '');
    outFile  = fullfile(outDir, sprintf('%s_panel%d.tex', scriptName, p));
    fid = fopen(outFile, 'w');
    fprintf(fid, '%s\n', panelLines{:});
    fclose(fid);
end
end
