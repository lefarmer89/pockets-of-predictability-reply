function results = run_smoke_tests(varargin)
% RUN_SMOKE_TESTS  Drive the package smoke-test suite.
%
%   results = run_smoke_tests()                  runs all tests
%   results = run_smoke_tests('only', {'sub'})   runs only matching tests
%   results = run_smoke_tests('verbose', true)   prints test stdout
%
% Each test is a function in this folder named test_*.m that returns a struct
% with fields:
%   .name     short test name
%   .pass     logical
%   .elapsed  seconds
%   .message  pass/fail summary string (one line)
%
% The driver prints a summary table and returns a struct array. Exit code
% (when run via -batch) is non-zero if any test fails.
%
% Convention: smoke tests must NOT trigger full production work. Subset data,
% cap iteration counts, or load cached .mat. Target total runtime under 5 min.

p = inputParser;
addParameter(p, 'only', {});
addParameter(p, 'verbose', false);
parse(p, varargin{:});
opts = p.Results;

here = fileparts(mfilename('fullpath'));
pkgRoot = fileparts(here);
addpath(pkgRoot);
% setup_paths registers the topical subroutine subdirs in one place
setup_paths();

testFiles = dir(fullfile(here, 'test_*.m'));
testNames = {testFiles.name};
testNames = cellfun(@(s) s(1:end-2), testNames, 'UniformOutput', false);

if ~isempty(opts.only)
    keep = false(size(testNames));
    for k = 1:numel(testNames)
        for j = 1:numel(opts.only)
            if contains(testNames{k}, opts.only{j})
                keep(k) = true; break
            end
        end
    end
    testNames = testNames(keep);
end

fprintf('\n=== Smoke tests ===\n');
fprintf('Package root: %s\n', pkgRoot);
fprintf('Running %d test(s)\n\n', numel(testNames));

results = repmat(struct('name','','pass',false,'elapsed',0,'message',''), ...
    1, numel(testNames));
for k = 1:numel(testNames)
    name = testNames{k};
    fprintf('  [%d/%d] %s ... ', k, numel(testNames), name);
    t0 = tic;
    try
        if opts.verbose
            r = feval(name);
        else
            % evalc captures stdout in arg 1; remaining args are the
            % expression's outputs. We discard stdout, keep r.
            [~, r] = evalc(sprintf('%s()', name));
        end
        r.elapsed = toc(t0);
        if ~isfield(r,'name'); r.name = name; end
        if r.pass
            fprintf('PASS (%.2fs)  %s\n', r.elapsed, r.message);
        else
            fprintf('FAIL (%.2fs)  %s\n', r.elapsed, r.message);
        end
    catch ME
        r = struct('name', name, 'pass', false, 'elapsed', toc(t0), ...
            'message', sprintf('ERROR %s: %s', ME.identifier, ME.message));
        fprintf('ERROR (%.2fs)  %s\n', r.elapsed, r.message);
    end
    results(k) = r;
end

passed = sum([results.pass]);
total  = numel(results);
totalTime = sum([results.elapsed]);
fprintf('\n=== Summary ===\n');
fprintf('  %d/%d passed in %.1fs\n', passed, total, totalTime);

if passed < total
    fprintf('  Failed tests:\n');
    for k = 1:total
        if ~results(k).pass
            fprintf('    - %s: %s\n', results(k).name, results(k).message);
        end
    end
    if isdeployed || ~usejava('desktop')
        exit(1);
    end
end
end
