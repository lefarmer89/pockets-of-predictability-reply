function computeTab1Results(cfg)
% COMPUTETAB1RESULTS  Aggregate the Fixed (sign-restricted, signSpec=2)
% and Adaptive (per-day expanding-window best by RMSE / by alpha)
% out-of-sample statistics into a single .mat file consumed by Table 1.
%
%   computeTab1Results(cfg)
%
% Reads:
%   results/daily/OOSResults_2_2.5yE_1yDM.mat
%       (statMat, statDiffMat, econMat — Fixed sign-restricted spec)
%   results/hyperparameters/forecastResults_2_2.5yE_1yDM_1S_pm_HyperR1.mat
%   results/hyperparameters/OOSResults_2_HyperR1_ExpandingC.mat
%       (consumed indirectly via adaptivePanelsForSplit with sampleSplit='full')
%
% Writes:
%   results/aggregates/tab1Results.mat with fields:
%     cwMatFixed         9 x 3   (FS, IP, OOP CW for Fixed)
%     econMatFixed       9 x 3   (alpha, t(alpha), SR for Fixed)
%     cwMatAdaptive      9 x 3 x 2  (CW for Adaptive {RMSE, alpha})
%     alphMatAdaptive    2 x 9   (annualized alpha by criterion)
%     tStatAlphMatAdaptive 2 x 9 (HAC t(alpha) by criterion)
%     SRAdaptive         2 x 9   (annualized SR by criterion)
%
% Used by: tables/table1.m
% Note: must run AFTER dailyEmpirics(cfg) and dailyEmpiricsHyperparameters{1,2}(cfg).

if nargin < 1 || isempty(cfg); cfg = default_config(); end
paths   = cfg.paths;
aggDir  = fullfile(paths.results, 'aggregates');
if ~isfolder(aggDir); mkdir(aggDir); end

%% Fixed (signSpec = 2, the published sign-restricted spec)
fixedFile = fullfile(paths.results, 'daily', 'OOSResults_2_2.5yE_1yDM.mat');
F = load(fixedFile, 'statMat', 'econMat');

% statMat order in dailyEmpirics: [FS, IP, OOP, Diff].
% Pull just the first 3 columns for Panel A; econMat already 3-col.
cwMatFixed   = F.statMat(:, 1:3);     %#ok<NASGU>
econMatFixed = F.econMat;             %#ok<NASGU>

%% Adaptive (per-day expanding-window best by RMSE and by alpha)
% adaptivePanelsForSplit returns:
%   cwBestMat[ii, c, a]      9 x 3 x 2
%   alphaAnnualBest[a, ii]   2 x 9
%   tStatAlphMatBest[a, ii]  2 x 9
%   SRBestAnnual[a, ii]      2 x 9
cfgAdaptive = cfg;
cfgAdaptive.sampleSplit = 'full';
[cwMatAdaptive, alphaAnnualBest, tStatAlphMatBest, SRBestAnnual] = ...
    adaptivePanelsForSplit(cfgAdaptive); %#ok<ASGLU>

% Match the variable names the legacy table1.m / Archive expect.
alphMatAdaptive      = alphaAnnualBest;     %#ok<NASGU>
tStatAlphMatAdaptive = tStatAlphMatBest;    %#ok<NASGU>
SRAdaptive           = SRBestAnnual;        %#ok<NASGU>

outFile = fullfile(aggDir, 'tab1Results.mat');
save(outFile, 'cwMatFixed', 'econMatFixed', 'cwMatAdaptive', ...
    'alphMatAdaptive', 'tStatAlphMatAdaptive', 'SRAdaptive');
fprintf('  wrote %s\n', outFile);
end
