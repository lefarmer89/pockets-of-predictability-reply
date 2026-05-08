function folder = robustnessResultsFolder(paths, cfg) %#ok<INUSD>
% ROBUSTNESSRESULTSFOLDER  Resolve the subfolder under paths.results
% where Hyperparameters{1,2,Marginals} cached outputs land.
%
%   folder = robustnessResultsFolder(paths, cfg)
%
% Returns paths.results/hyperparameters. cfg is accepted for call-site
% symmetry with baselineResultsFolder but is not currently consumed
% (sample splits are now handled ex-post via computeAdaptiveSampleSplit;
% no fresh sub-sample sweeps are written).
%
% The folder is created if it does not exist.

folder = fullfile(paths.results, 'hyperparameters');

if ~isfolder(folder)
    mkdir(folder);
end
end
