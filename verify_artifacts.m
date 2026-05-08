function passed = verify_artifacts()
%VERIFY_ARTIFACTS  Walk results/ and confirm every cached .mat referenced
% by the reply tables and figures is present.
%
%   passed = verify_artifacts()
%
% Prints a PASS/FAIL line per expected artifact and returns true iff every
% expected file exists.

setup_paths();
paths = struct('results', fullfile(pwd, 'results'));

% Daily empirics: per-spec forecastResults / OOSResults plus aggregates.
checks = struct('relPath', {}, 'requiredVars', {});

dailySpecs = {
    {0, 0, 2.7, 1};   {0, 0, 2.5, 0.5}; {0, 0, 2.5, 1.5};
    {0, 0, 2,   1};   {0, 0, 3,   1};
    {1, 0, 2.7, 1};   {1, 0, 2.5, 0.5}; {1, 0, 2.5, 1.5};
    {1, 0, 2,   1};   {1, 0, 3,   1};
    {1, 1, 2.7, 1};   {1, 1, 2.5, 0.5}; {1, 1, 2.5, 1.5};
    {1, 1, 2,   1};   {1, 1, 3,   1}};
for k = 1:numel(dailySpecs)
    spec = dailySpecs{k};
    [eLabel, dmLabel, rLabel] = buildSpecLabels(cell2mat(spec));
    forecastFile = sprintf('forecastResults%s%s%s_1S_pm.mat', rLabel, eLabel, dmLabel);
    oosFile      = sprintf('OOSResults%s%s%s.mat',       rLabel, eLabel, dmLabel);
    checks(end+1) = struct('relPath', fullfile('daily', forecastFile), ...
        'requiredVars', {{'yF1Mat', 'yF2Mat', 'pocketIndMat', 'integralR2Mat', 'durationMat'}}); %#ok<AGROW>
    checks(end+1) = struct('relPath', fullfile('daily', oosFile), ...
        'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}}); %#ok<AGROW>
end

% Monthly empirics
monthlySpecs = {{0, 0}, {1, 0}, {1, 1}};
for k = 1:numel(monthlySpecs)
    sr = monthlySpecs{k}{1}; cr = monthlySpecs{k}{2};
    sgn = sr*1 + cr*1 + 1;   % 1, 2, 3
    checks(end+1) = struct('relPath', fullfile('monthly', sprintf( ...
        'forecastResultsMonthly_%d_2.5yE_1yDM_1S_pm.mat', sgn)), ...
        'requiredVars', {{'durationMat', 'integralR2Mat', 'pocketIndMat'}}); %#ok<AGROW>
    checks(end+1) = struct('relPath', fullfile('monthly', sprintf( ...
        'OOSResultsMonthly_%d_2.5yE_1yDM.mat', sgn)), ...
        'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}}); %#ok<AGROW>
end

% Fama-French
checks(end+1) = struct('relPath', fullfile('famaFrench', 'OOSResultsSMB_1_2.5yE_1yDM.mat'), ...
    'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}});
checks(end+1) = struct('relPath', fullfile('famaFrench', 'OOSResultsHML_1_2.5yE_1yDM.mat'), ...
    'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}});

% G-prior robustness
for g = [1, 3]
    checks(end+1) = struct('relPath', fullfile('gPriorRobustness', sprintf('g%d', g), ...
        'OOSResults_2_2.5yE_1yDM.mat'), ...
        'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}}); %#ok<AGROW>
end

% Hyperparameter sweep (signSpec ∈ {1, 2})
for s = [1, 2]
    checks(end+1) = struct('relPath', fullfile('hyperparameters', sprintf( ...
        'forecastResults_%d_2.5yE_1yDM_1S_pm_HyperR1.mat', s)), ...
        'requiredVars', {{'yActual', 'yF1PMMat', 'yF2All', 'pocketIndAll', 'paramCombs'}}); %#ok<AGROW>
    checks(end+1) = struct('relPath', fullfile('hyperparameters', sprintf( ...
        'OOSResults_%d_HyperR1_ExpandingC.mat', s)), ...
        'requiredVars', {{'alphMatExpanding', 'tStatAlphMatExpanding'}}); %#ok<AGROW>
end

% Aggregates
checks(end+1) = struct('relPath', fullfile('aggregates', 'tab1Results.mat'), ...
    'requiredVars', {{'cwMatFixed', 'cwMatAdaptive', 'econMatFixed', ...
                      'alphMatAdaptive', 'tStatAlphMatAdaptive', 'SRAdaptive'}});
checks(end+1) = struct('relPath', fullfile('aggregates', 'tabA3Results.mat'), ...
    'requiredVars', {{'econMatFixed', 'econMatFixed5', 'econMatFixed10', ...
                      'alphMatAdaptive', 'alphMatAdaptive5', 'alphMatAdaptive10'}});
checks(end+1) = struct('relPath', fullfile('aggregates', 'topKbotK.mat'), ...
    'requiredVars', {{'tStatAlphMatTopK', 'cwInMatTopK', 'tStatAlphaBenchmark', 'cwBenchmark'}});
checks(end+1) = struct('relPath', fullfile('aggregates', 'signifCell.mat'), ...
    'requiredVars', {{'signifCell'}});

% Asset-pricing model bootstraps
apModels = {'BY', 'CC', 'GP', 'W', 'W_nd'};
for k = 1:numel(apModels)
    checks(end+1) = struct('relPath', fullfile('assetPricing', sprintf( ...
        '%s_asset_pricing_sims_signRestriction_0_coefRestriction_0_25ybandwidth_1yDM_OOS.mat', apModels{k})), ...
        'requiredVars', {{'cwMat', 'cwDiffMat', 'econMat', 'tStatAlphaMat'}}); %#ok<AGROW>
end

% Sticky-expectations sims
for prefix = {'SE', 'RE', 'RE_Recalibrated'}
    checks(end+1) = struct('relPath', fullfile('assetPricing', sprintf( ...
        '%s_asset_pricing_sims_signRestriction_0_coefRestriction_0_25ybandwidth_1yDM_OOS.mat', prefix{1})), ...
        'requiredVars', {{'cwMat', 'cwDiffMat', 'econMat', 'tStatAlphaMat'}}); %#ok<AGROW>
end

% Run the checks
nPass = 0;
nFail = 0;
for k = 1:numel(checks)
    relPath = checks(k).relPath;
    absPath = fullfile(paths.results, relPath);
    if ~isfile(absPath)
        fprintf('FAIL  missing file: %s\n', relPath);
        nFail = nFail + 1;
        continue
    end
    info = whos('-file', absPath);
    presentVars = {info.name};
    missingVars = setdiff(checks(k).requiredVars, presentVars);
    if ~isempty(missingVars)
        fprintf('FAIL  %s: missing var(s) %s\n', relPath, strjoin(missingVars, ', '));
        nFail = nFail + 1;
        continue
    end
    fprintf('PASS  %s\n', relPath);
    nPass = nPass + 1;
end

fprintf('\n--- verify_artifacts: %d PASS, %d FAIL ---\n', nPass, nFail);
passed = (nFail == 0);
end
