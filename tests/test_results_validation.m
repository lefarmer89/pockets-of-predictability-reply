function r = test_results_validation()
% TEST_RESULTS_VALIDATION  Confirm that pre-computed .mat result files exist
% and have the expected variable structure.
%
% This guards against accidental deletion or corruption of the cached
% pipeline outputs that the table*.m and figure*.m scripts consume. Only
% reads .mat files; does not run any compute.

r.name = 'test_results_validation';
r.pass = false;
r.message = '';

here = fileparts(mfilename('fullpath'));
pkgRoot = fileparts(here);

required = {
    {'results/aggregates/tab1Results.mat',  {'cwMatFixed','cwMatAdaptive','econMatFixed', ...
                                             'alphMatAdaptive','tStatAlphMatAdaptive','SRAdaptive'}};
    {'results/aggregates/tabA3Results.mat', {'econMatFixed','alphMatAdaptive','tStatAlphMatAdaptive', ...
                                             'SRAdaptive','econMatFixed5','econMatFixed10'}};
    {'data/pcPredictor.mat',                {'predictors'}};
    {'results/aggregates/signifCell.mat',   {'signifCell'}};
    {'results/daily/forecastResults_1_2.5yE_1yDM_1S_pm.mat', ...
                                            {'durationMat','integralR2Mat','pocketIndMat'}};
};

found = 0;
missing = {};
varMissing = {};

for k = 1:numel(required)
    fname = required{k}{1};
    expectedVars = required{k}{2};
    fpath = fullfile(pkgRoot, fname);
    if ~isfile(fpath)
        missing{end+1} = fname; %#ok<AGROW>
        continue
    end
    found = found + 1;
    try
        info = whos('-file', fpath);
        actualVars = {info.name};
        for j = 1:numel(expectedVars)
            if ~ismember(expectedVars{j}, actualVars)
                varMissing{end+1} = sprintf('%s::%s', fname, expectedVars{j}); %#ok<AGROW>
            end
        end
    catch ME
        varMissing{end+1} = sprintf('%s (read error: %s)', fname, ME.message); %#ok<AGROW>
    end
end

if isempty(missing) && isempty(varMissing)
    r.pass = true;
    r.message = sprintf('all %d required artifacts present with expected variables', numel(required));
else
    msgs = {};
    if ~isempty(missing)
        msgs{end+1} = sprintf('missing files: %s', strjoin(missing, ', '));
    end
    if ~isempty(varMissing)
        msgs{end+1} = sprintf('missing vars: %s', strjoin(varMissing, ', '));
    end
    r.message = strjoin(msgs, '; ');
end
end
