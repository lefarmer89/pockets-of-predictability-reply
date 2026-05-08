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

% Daily empirics: canonical forecastResults plus the 7 OOS files that
% displayResults consumes (table3 needs sgn=1/2/3 canonical, tableA2 needs
% the sgn=2 window-robustness sweep, table2/figure2/computeSignifCell need
% the canonical sgn=1 forecast).
checks = struct('relPath', {}, 'requiredVars', {});

checks(end+1) = struct('relPath', fullfile('daily', 'forecastResults_1_2.5yE_1yDM_1S_pm.mat'), ...
    'requiredVars', {{'yF1Mat', 'yF2Mat', 'pocketIndMat', 'integralR2Mat', 'durationMat'}});

dailyOosFiles = {
    'OOSResults_1_2.5yE_1yDM.mat'    % sgn=1 canonical (table3)
    'OOSResults_2_2.5yE_1yDM.mat'    % sgn=2 canonical (table3, aggregates)
    'OOSResults_2_2.5yE_6mDM.mat'    % sgn=2 windowed (tableA2)
    'OOSResults_2_2.5yE_1.5yDM.mat'  % sgn=2 windowed (tableA2)
    'OOSResults_2_2yE_1yDM.mat'      % sgn=2 windowed (tableA2)
    'OOSResults_2_3yE_1yDM.mat'      % sgn=2 windowed (tableA2)
    'OOSResults_3_2.5yE_1yDM.mat'};  % sgn=3 canonical (table3)
for k = 1:numel(dailyOosFiles)
    checks(end+1) = struct('relPath', fullfile('daily', dailyOosFiles{k}), ...
        'requiredVars', {{'statMat', 'statDiffMat', 'econMat'}}); %#ok<AGROW>
end

% Monthly empirics: only the sgn=1 forecast file (consumed by tableA5);
% all 3 OOS files (consumed by tableA6).
checks(end+1) = struct('relPath', fullfile('monthly', ...
    'forecastResultsMonthly_1_2.5yE_1yDM_1S_pm.mat'), ...
    'requiredVars', {{'durationMat', 'integralR2Mat', 'pocketIndMat'}});
for sgn = 1:3
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

% Hyperparameter sweep intermediates (forecastResults_*_HyperR1.mat,
% OOSResults_*_HyperR1_ExpandingC.mat) are produced by Hyperparameters1/2,
% consumed by HyperparametersMarginals to build aggregates/topKbotK.mat,
% and then no longer needed. They are intentionally absent from the public
% archive; verify_artifacts checks the downstream aggregates only.

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

% Asset-pricing model bootstraps (table4 consumes BY/CC/DT/GP/W/W_nd)
apModels = {'BY', 'CC', 'DT', 'GP', 'W', 'W_nd'};
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
