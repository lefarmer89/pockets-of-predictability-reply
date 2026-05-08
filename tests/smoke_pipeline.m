function smoke_pipeline()
%SMOKE_PIPELINE  End-to-end smoke test of the production pipeline.
%
%   smoke_pipeline()
%
% Calls every production entry point with shakeout configuration and
% redirects all .mat writes to results_smoke/ so existing cached
% results/ files are NOT overwritten.
%
% Shakeout config:
%   cfg.dailySpecsSubset           = 1   (1 of 15 specs in dailyEmpirics)
%   cfg.monthlySpecsSubset         = 1   (1 of 3 specs in monthlyEmpirics)
%   cfg.hyperparameterComboSubset  = 1:10  (10 of 9,720 combos)
%   cfg.hyperparameterSpecs        = 2   (1 of 2 signSpecs)
%   cfg.bootstrapReps              = 5
%
% Reports per-stage timing and a final PASS / FAIL line per stage. A
% stage that errors continues so we can see which stages survive in
% isolation.

setup_paths();
cfg = default_config();

% Redirect all cached outputs to a scratch directory.
cfg.paths.results = fullfile(pwd, 'results_smoke');
if ~isfolder(cfg.paths.results); mkdir(cfg.paths.results); end

% Shakeout knobs.
cfg.dailySpecsSubset          = 1;
cfg.monthlySpecsSubset        = 1;
cfg.hyperparameterComboSubset = 1:10;
cfg.hyperparameterSpecs       = 2;
cfg.bootstrapReps             = 5;

% Stage list: {label, callback}. Each callback takes cfg.
stages = {
    'dailyEmpirics',                    @(cfg) dailyEmpirics(cfg);
    'dailyEmpiricsGPriorRobustness',    @(cfg) dailyEmpiricsGPriorRobustness(cfg);
    'monthlyEmpirics',                  @(cfg) monthlyEmpirics(cfg);
    'famaFrenchEmpirics',               @(cfg) famaFrenchEmpirics(cfg, {'SMB'});
    'dailyEmpiricsHyperparameters1',    @(cfg) dailyEmpiricsHyperparameters1(cfg);
    'dailyEmpiricsHyperparameters2',    @(cfg) dailyEmpiricsHyperparameters2(cfg, 2);
    'dailyEmpiricsHyperparametersMarginals', @(cfg) dailyEmpiricsHyperparametersMarginals(cfg, 2);
    'dailyBootstrap',                   @(cfg) dailyBootstrap(cfg);
    'generateStickyExpectationsData',   @(cfg) generateStickyExpectationsData(cfg);
    'stickyExpectationsSim',            @(cfg) stickyExpectationsSim(cfg);
    'dailyAssetPricingBootstrap',       @(cfg) dailyAssetPricingBootstrap(cfg);
    'computeSignifCell',                @(cfg) computeSignifCell(cfg);
};
% computeTab1Results and computeTabA3Results are aggregate consumers
% that read multi-spec daily/ outputs (specs 1, 6, 11). With
% cfg.dailySpecsSubset = 1, those files are not produced; skip rather
% than chase aggregator-only failures unrelated to the pipeline.

results = struct('label', {}, 'elapsed', {}, 'status', {}, 'message', {});
totalStart = tic;

fprintf('\n=== smoke_pipeline ===\n');
fprintf('Output dir: %s\n', cfg.paths.results);
fprintf('Stages:     %d\n\n', size(stages, 1));

for k = 1:size(stages, 1)
    label = stages{k, 1};
    cb    = stages{k, 2};
    fprintf('[%2d/%2d] %s ... ', k, size(stages, 1), label);
    tStart = tic;
    try
        cb(cfg);
        elapsed = toc(tStart);
        fprintf('PASS (%6.1fs)\n', elapsed);
        results(end+1) = struct('label', label, 'elapsed', elapsed, ...
            'status', 'PASS', 'message', ''); %#ok<AGROW>
    catch ME
        elapsed = toc(tStart);
        fprintf('FAIL (%6.1fs): %s\n', elapsed, ME.message);
        results(end+1) = struct('label', label, 'elapsed', elapsed, ...
            'status', 'FAIL', 'message', ME.message); %#ok<AGROW>
    end
end

totalElapsed = toc(totalStart);

fprintf('\n=== Summary ===\n');
nPass = sum(strcmp({results.status}, 'PASS'));
nFail = sum(strcmp({results.status}, 'FAIL'));
fprintf('  %d / %d PASS in %.1f min\n', nPass, nFail + nPass, totalElapsed / 60);
for k = 1:numel(results)
    fprintf('    %-42s %s (%6.1fs)\n', results(k).label, results(k).status, results(k).elapsed);
end
end
