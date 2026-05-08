function s = prettifyLatexTable(latexStr)
%PRETTIFYLATEXTABLE  Convert LatexTableFull output into a readable form.
%
%   s = prettifyLatexTable(latexStr)
%
% Strips LaTeX preamble/postamble lines and inline markup ($, \cr, \\,
% \textbf{}, \hat{}, ^{*}, ^{**}, ^{***}, ^{\dagger}, etc.), and replaces
% column separators (` & `) with vertical bars so the table renders
% cleanly in the MATLAB console.
%
% Designed for the output of subroutines/io/LatexTableFull.m.

lines = regexp(latexStr, '\r\n|\n|\r', 'split');
out = cell(1, numel(lines));
for k = 1:numel(lines)
    line = lines{k};

    % Drop LaTeX preamble/postamble lines.
    if startsWith(strtrim(line), '\') && ...
            (contains(line, {'\newpage', '\begin', '\end', '\centering', ...
                             '\caption', '\resizebox', '\\'}))
        out{k} = '';
        continue
    end

    % Strip end-of-row \cr.
    line = regexprep(line, '\s*\\cr\s*$', '');

    % Strip $...$ math wrappers.
    line = regexprep(line, '\$([^$]*)\$', '$1');

    % Stars and daggers (significance markers).
    line = regexprep(line, '\^\{\\dagger\\dagger\\dagger\}', '+++');
    line = regexprep(line, '\^\{\\dagger\\dagger\}', '++');
    line = regexprep(line, '\^\{\\dagger\}', '+');
    line = regexprep(line, '\^\{\*\*\*\}', '***');
    line = regexprep(line, '\^\{\*\*\}', '**');
    line = regexprep(line, '\^\{\*\}', '*');
    line = regexprep(line, '\^\{[^}]*\}', '');

    % Subscripts: t_{\alpha} → t_alpha (drop braces, keep contents).
    line = regexprep(line, '_\{([^}]*)\}', '_$1');

    % Inline LaTeX commands (drop wrapper, keep argument).
    line = regexprep(line, '\\hat\{([^}]*)\}', '$1');
    line = regexprep(line, '\\textbf\{([^}]*)\}', '$1');
    line = regexprep(line, '\\multicolumn\{[^}]*\}\{[^}]*\}\{([^}]*)\}', '$1');
    line = regexprep(line, '\\[a-zA-Z]+\{([^}]*)\}', '$1');

    % Bare LaTeX command names: \alpha → alpha, \beta → beta, etc.
    line = regexprep(line, '\\([a-zA-Z]+)', '$1');

    % Drop any leftover braces.
    line = strrep(line, '{', '');
    line = strrep(line, '}', '');

    % Replace column separators with vertical bars. Also strip a trailing
    % `&` (no whitespace required after) so end-of-row reads cleanly.
    line = regexprep(line, '\s+&\s*$', '');
    line = regexprep(line, '\s+&\s+', '  | ');

    out{k} = line;
end

% Drop empty strings created by preamble removal.
out = out(~cellfun('isempty', strtrim(out)));
s = strjoin(out, newline);
end
