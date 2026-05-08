function computeTabA3Results(cfg)
% COMPUTETABA3RESULTS  Aggregate the Fixed (sign-restricted, signSpec=2)
% and Adaptive (per-day expanding-window best by RMSE / by alpha)
% economic-significance statistics at three transaction-cost levels
% (0, 5, 10 basis points) into a single .mat file consumed by Table A.3.
%
%   computeTabA3Results(cfg)
%
% Reads (assumes runAll has already produced these for each bps level):
%   results/daily_<bps>bps/OOSResults_2_2.5yE_1yDM.mat
%       (econMat at adjCostBps = bps)        — Fixed columns
%   results/aggregates/topKbotK_<bps>bps.mat
%       (alphMatBest, tStatAlphMatBest, SRBest at adjCostBps = bps) — Adaptive columns
% Falls back to results/daily/ and results/aggregates/topKbotK.mat for
% the 0 bps slot to remain consistent with the published baseline.
%
% Writes:
%   results/aggregates/tabA3Results.mat with fields:
%     econMatFixed, alphMatAdaptive, tStatAlphMatAdaptive, SRAdaptive       (0 bps)
%     econMatFixed5, alphMatAdaptive5, tStatAlphMatAdaptive5, SRAdaptive5   (5 bps)
%     econMatFixed10, alphMatAdaptive10, tStatAlphMatAdaptive10, SRAdaptive10 (10 bps)
%
% Used by: tables/tableA3.m
% Note: must run AFTER 3 successive runAll(cfg) executions with
% cfg.adjCostBps set to 0, 5, 10 respectively. See runAll.m for the
% staging helper that loops over bps levels and writes the bps-suffixed
% artifacts that this function consumes.

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths   = cfg.paths;
aggDir  = fullfile(paths.results, 'aggregates');
if ~isfolder(aggDir); mkdir(aggDir); end

bpsLevels = [0, 5, 10];
suffixes  = {'', '5', '10'};

allOut = struct();
for ll = 1:numel(bpsLevels)
    bps = bpsLevels(ll);
    [econMatFixed, alphMatAdaptive, tStatAlphMatAdaptive, SRAdaptive] = ...
        oneBpsLevel(cfg, bps);

    s = suffixes{ll};
    allOut.(['econMatFixed' s])         = econMatFixed;
    allOut.(['alphMatAdaptive' s])      = alphMatAdaptive;
    allOut.(['tStatAlphMatAdaptive' s]) = tStatAlphMatAdaptive;
    allOut.(['SRAdaptive' s])           = SRAdaptive;
end

outFile = fullfile(aggDir, 'tabA3Results.mat');
save(outFile, '-struct', 'allOut');
fprintf('  wrote %s\n', outFile);
end


function [econMatFixed, alphMatAdaptive, tStatAlphMatAdaptive, SRAdaptive] = ...
    oneBpsLevel(cfg, bps)
% Pull Fixed and Adaptive results at the specified transaction-cost
% level. For bps = 0 we read from the canonical (un-suffixed) folders;
% for 5 / 10 we read from the bps-suffixed ones produced by the
% transaction-cost loop in runAll.

paths = cfg.paths;
if bps == 0
    dailySub = 'daily';
    aggSfx   = '';
else
    dailySub = sprintf('daily_%dbps', bps);
    aggSfx   = sprintf('_%dbps', bps);
end

%% Fixed
fixedFile = fullfile(paths.results, dailySub, 'OOSResults_2_2.5yE_1yDM.mat');
F = load(fixedFile, 'econMat');
econMatFixed = F.econMat;

%% Adaptive
% adaptivePanelsForSplit caches its own input; rerun it with the
% transaction-cost-aware cfg so the per-period regressions are evaluated
% at the bps-adjusted yPocketTime. NOTE: this requires the H1/H2 caches
% themselves to have been regenerated at the bps level, which the
% transaction-cost stage of runAll arranges via dailyEmpiricsHyperparameters{1,2}.
cfgAdaptive = cfg;
cfgAdaptive.sampleSplit = 'full';
cfgAdaptive.adjCostBps  = bps;
% Look for a bps-suffixed cache first; else fall back to recompute.
cwBestCache = fullfile(paths.results, 'aggregates', sprintf('cwBestMat_full%s.mat', aggSfx));
if isfile(cwBestCache)
    C = load(cwBestCache);
    alphMatAdaptive      = C.alphMatAdaptive;
    tStatAlphMatAdaptive = C.tStatAlphMatAdaptive;
    SRAdaptive           = C.SRAdaptive;
else
    [~, alphMatAdaptive, tStatAlphMatAdaptive, SRAdaptive] = ...
        adaptivePanelsForSplit(cfgAdaptive);
end
end
