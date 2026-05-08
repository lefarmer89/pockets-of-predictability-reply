function r = test_yF2AllS_bugfix()
% TEST_YF2ALLS_BUGFIX  Verify that the yF2AllS construction in the
% refactored dailyEmpiricsHyperparameters2 reproduces the archive's
% yF2AllS values within 1e-10. The previous refactor wrongly computed
% yComb2 / yComb3 from the in-pocket-mutated yF2All; the archive
% computes them from the ORIGINAL yF2All (BEFORE mutation), and yComb1
% from the MUTATED yF2All (AFTER mutation).
%
% Test strategy: load the inputs (yF2All, pocketIndAll, yF1PMMat) from
% the current forecastResults file, run them through the new yF2AllS
% construction, and compare to the archive's yF2AllS for the first 10
% combos (a small enough slice to load quickly from a v7.3 .mat).

r = struct('name', mfilename, 'pass', false, 'elapsed', 0, 'message', '');

here     = fileparts(mfilename('fullpath'));
pkgRoot  = fileparts(here);
projRoot = fileparts(pkgRoot);
addpath(pkgRoot);

curFile = fullfile(pkgRoot, 'results', 'hyperparameters', ...
    'forecastResults_2_2.5yE_1yDM_1S_pm_HyperR1.mat');
arcFile = fullfile(projRoot, 'Archive', 'Replication Package - Original Backup', ...
    'Robustness_Results', 'forecastResults_2_2.5yE_1yDM_1S_pm_HyperR1.mat');

if ~isfile(curFile);  r.message = sprintf('current file missing: %s', curFile); return; end
if ~isfile(arcFile);  r.message = sprintf('archive file missing: %s', arcFile); return; end

% Slice the inputs to the first 10 combos via matfile (no full load).
mC = matfile(curFile);
yF2AllSlice       = mC.yF2All(:, :, 1:10) ./ 100;       % T_full x 6 x 10
pocketIndAllSlice = mC.pocketIndAll(:, :, 1:10);         % T_full x 6 x 10
yF1PMMat          = mC.yF1PMMat ./ 100;                  % T_full x 6
yActual           = mC.yActual ./ 100;                   % T_full x 1

% Trim. The trim value is computed across the FULL grid in production
% (max over all 9720 combos); subset would shrink it. To match the
% archive's row alignment, recompute trim from the full file.
trim = max(sum(isnan(mC.yF2All)), [], 'all');
yF1PMMat        = yF1PMMat(trim+1:end, :);
yF2AllSlice     = yF2AllSlice(trim+1:end, :, :);
yActual         = yActual(trim+1:end, :);
pocketIndAllSlice = pocketIndAllSlice(trim+1:end, :, :);
pocketIndAllSlice(isnan(pocketIndAllSlice)) = 0;
pocketIndAllSlice = logical(pocketIndAllSlice);

% Compute yF2AllS the NEW (bugfixed) way.
yF2All       = yF2AllSlice;       % alias for clarity
pocketIndAll = pocketIndAllSlice;
yComb2 = combineInPocketAverageLocal(yF2All, pocketIndAll, yF1PMMat);
yComb3 = squeeze(mean(yF2All(:, 1:4, :), 2, 'omitmissing'));
yF2InPocket = yF2All .* pocketIndAll + yF1PMMat .* (~pocketIndAll);
yComb1 = squeeze(mean(yF2InPocket(:, 1:4, :), 2, 'omitmissing'));

T = size(yActual, 1);
yF2AllS_new = NaN(T, size(yF2All, 2) + 3, size(yF2All, 3));
yF2AllS_new(:, 1:6, :) = yF2InPocket;
yF2AllS_new(:, 7, :)   = yComb1;
yF2AllS_new(:, 8, :)   = yComb2;
yF2AllS_new(:, 9, :)   = yComb3;

% Load archive's yF2AllS slice for combos 1..10. (Archive stores yF2AllS
% in this same forecastResults file, appended by the original H2.)
mA = matfile(arcFile);
yF2AllS_archive = mA.yF2AllS(:, :, 1:10);

% Compare.
diff = abs(yF2AllS_new - yF2AllS_archive);
maxDiff = max(diff(:), [], 'omitmissing');
nanMismatch = sum(isnan(yF2AllS_new(:)) ~= isnan(yF2AllS_archive(:)));

% Per-column diagnostic (verbose).
for col = 1:9
    sub = yF2AllS_new(:, col, :) - yF2AllS_archive(:, col, :);
    md  = max(abs(sub(:)), [], 'omitmissing');
    fprintf('  col %d:  maxAbsDiff = %.3g\n', col, md);
end

r.pass    = (maxDiff < 1e-10) && (nanMismatch == 0);
r.message = sprintf('maxAbsDiff = %.3g  nanMismatch = %d', maxDiff, nanMismatch);
end


function yComb2 = combineInPocketAverageLocal(yF2, pocketIndAll, yF1PMMat)
% Local copy of dailyEmpiricsHyperparameters2's combineInPocketAverage.
T = size(yF2, 1);
numParamCombs = size(yF2, 3);
yComb2 = NaN(T, numParamCombs);
for ii = 1:T
    for jj = 1:numParamCombs
        if any(pocketIndAll(ii, 1:4, jj))
            yComb2(ii, jj) = mean(yF2(ii, logical(pocketIndAll(ii, 1:4, jj)), jj));
        else
            yComb2(ii, jj) = yF1PMMat(ii, 1);
        end
    end
end
end
