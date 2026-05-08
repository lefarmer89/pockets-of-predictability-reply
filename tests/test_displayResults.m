function r = test_displayResults()
% TEST_DISPLAYRESULTS  Run displayResults end-to-end on cached .mat files
% and confirm it completes without error. Captures stdout to a log file.
%
% Only consumes existing pipeline outputs - does NOT trigger compute. The
% existing displayResults / table*.m / figure*.m use `clearvars`, which
% would wipe the local result struct - so we run the orchestrator in a
% subprocess and check exit code.

r.name = 'test_displayResults';
r.pass = false;
r.message = '';

here = fileparts(mfilename('fullpath'));
pkgRoot = fileparts(here);

% Find a MATLAB executable
matlabExe = fullfile(matlabroot, 'bin', 'matlab.exe');
if ~isfile(matlabExe)
    matlabExe = fullfile(matlabroot, 'bin', 'matlab');
    if ~isfile(matlabExe)
        r.message = 'cannot locate matlab executable';
        return
    end
end

logFile = fullfile(pkgRoot, 'output', 'logs', ...
    sprintf('test_displayResults_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
if ~isfolder(fileparts(logFile))
    mkdir(fileparts(logFile));
end

% Build command: cd into pkgRoot, run displayResults, close all, exit
cmd = sprintf('"%s" -batch "cd(''%s''); displayResults; close all; exit"', ...
    matlabExe, strrep(pkgRoot, '\', '\\'));
cmd = sprintf('%s 1> "%s" 2>&1', cmd, logFile);

t0 = tic;
[status, ~] = system(cmd);
elapsed = toc(t0);

if status == 0
    r.pass = true;
    r.message = sprintf('displayResults completed in %.1fs (log: %s)', ...
        elapsed, logFile);
else
    % Read tail of log for diagnosis
    tail = '';
    if isfile(logFile)
        fid = fopen(logFile, 'r');
        if fid > 0
            txt = fread(fid, '*char')';
            fclose(fid);
            lines = strsplit(txt, newline);
            keep = max(1, numel(lines)-3);
            tail = strjoin(lines(keep:end), ' | ');
        end
    end
    r.message = sprintf('displayResults exit=%d, log=%s, tail: %s', ...
        status, logFile, tail);
end
end
