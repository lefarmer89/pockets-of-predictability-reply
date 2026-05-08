function generateStickyExpectationsData(cfg)
% GENERATESTICKYEXPECTATIONSDATA  Simulate three asset-pricing model paths
% (sticky, rational, rational-recalibrated) and write them to
% results/simulatedPaths/. Inputs to stickyExpectationsSim.m and to
% Section V of the reply.
%
% Inputs (read from cfg.paths.data):
%   stickyCalibration_freevols_r211.mat   (theta vector for sticky model)
%   rationalCalibration_freevols_r211.mat (theta vector for rational recalibration)
%
% Outputs (written to cfg.paths.results/simulatedPaths/):
%   sticky_sims.mat
%   rational_sims.mat
%   rational_recalibrated_sims.mat
%
% Tables/figures consumed by:
%   Table 5 (via stickyExpectationsSim.m)
%
% Runtime: ~10-15 minutes per simulation on a single core.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

paths = cfg.paths;

cal = ["sticky", "rational", "recalibrated"];
TSample = 12000000;

outDir = fullfile(paths.results, 'simulatedPaths');
if ~isfolder(outDir); mkdir(outDir); end

stickyCal = load(fullfile(paths.data, 'stickyCalibration_freevols_r211.mat'));
rationalCal = load(fullfile(paths.data, 'rationalCalibration_freevols_r211.mat'));

tic
for i = 1:3
    spec = cal(i);
    fprintf('Simulating %s\n', spec);

    switch spec
        case "sticky"
            theta = stickyCal.thetaStickyNew;
            lambda = 0.3^(4/252);
            outMatfile = fullfile(outDir, 'sticky_sims.mat');
        case "rational"
            theta = stickyCal.thetaStickyNew;
            lambda = 0;
            outMatfile = fullfile(outDir, 'rational_sims.mat');
        case "recalibrated"
            theta = rationalCal.thetaRationalNew;
            lambda = 0;
            outMatfile = fullfile(outDir, 'rational_recalibrated_sims.mat');
    end

    PathSimulator(theta, TSample, cfg.rngSeed, lambda, outMatfile);
end
toc
end
