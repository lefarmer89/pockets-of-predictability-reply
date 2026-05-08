function fp = saveTableScript(funcName, ~, outDir)
%SAVETABLESCRIPT  Run a table-rendering function with stdout capture, save
% the captured LaTeX to a .tex file, and print a stripped readable form.
%
%   fp = saveTableScript(funcName, [unused], outDir)
%
% funcName  e.g. 'table1', 'tableA3'. A function in tables/ that calls
%           LatexTableFull (which writes via disp).
% outDir    Output directory for the .tex file. Defaults to
%           'output/tables' relative to the package root.
%
% Returns the absolute path of the written .tex file.

if nargin < 3 || isempty(outDir); outDir = fullfile('output', 'tables'); end
if ~isfolder(outDir); mkdir(outDir); end

fp = fullfile(outDir, [funcName, '.tex']);
captured = evalc(funcName);

fid = fopen(fp, 'w');
fprintf(fid, '%s', captured);
fclose(fid);

% Emit a data-only .tex (no header rows) for use via \input{} from the
% reply LaTeX. Drop lines that begin with the column-header sentinels
% LatexTableFull emits via the rowString{1} slot, plus any blank lines.
% Rewrite the plain-TeX `\cr` row terminator that LatexTableFull uses to
% LaTeX `\\`. Strip the trailing `\\` from the last data line so the
% parent reply LaTeX can supply its own row terminator after \input{...}.
fpData = strrep(fp, '.tex', '_data.tex');
lines = regexp(captured, '\r?\n', 'split');
isHeader = @(L) startsWith(strtrim(L), 'Variables &') || ...
                startsWith(strtrim(L), 'Statistics &') || ...
                startsWith(strtrim(L), 'Stats &') || ...
                startsWith(strtrim(L), 'Data &');
mask  = cellfun(@(L) ~isempty(strtrim(L)) && ~isHeader(L), lines);
linesData = strrep(lines(mask), '\cr', '\\');
linesData{end} = regexprep(linesData{end}, '\s*\\\\\s*$', '');
fid = fopen(fpData, 'w');
fprintf(fid, '%s\n', linesData{:});
fclose(fid);

fprintf('\n=== %s  →  %s ===\n', funcName, fp);
fprintf('%s\n', prettifyLatexTable(captured));
end
