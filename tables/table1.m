function table1()
%TABLE1  Render Table 1: out-of-sample forecasting performance.
%
% Panel A reports Clark-West statistics, Panel B reports alpha, t(alpha),
% and Sharpe ratio. The five column groups are FST (2023), CFNPZ (2024),
% Fixed, Adaptive (min RMSE), and Adaptive (max alpha). FST and CFNPZ
% values are reproduced verbatim from the source manuscripts; the other
% three groups load from results/aggregates/tab1Results.mat.
%
% Reads:  results/aggregates/tab1Results.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

load('results/aggregates/tab1Results.mat', 'cwMatFixed', 'cwMatAdaptive', ...
    'econMatFixed', 'alphMatAdaptive', 'tStatAlphMatAdaptive', 'SRAdaptive');

% Rows: dp, tbl, tsp, rvar, mv, pc, comb1, comb2, comb3.
% Panel A columns: full sample CW, in-pocket CW, out-of-pocket CW.
fstA = [
     0.40,  3.79, -1.94;
     1.98,  4.75, -1.33;
     0.95,  4.52, -1.54;
    -0.79,  3.93, -1.07;
    -0.01,  4.01, -1.22;
     1.85,  4.69, -0.52;
     6.35,  6.55,  NaN;
     6.15,  6.33,  NaN;
     0.23,  0.83, -1.23];

cfnpzA = [
     0.40, -0.64,  0.87;
     1.98,  2.02,  0.52;
     0.95,  0.90,  0.44;
    -0.79,  0.38, -0.93;
    -0.01,  0.23, -0.29;
     1.85,  2.19,  0.52;
     1.12,  1.12,  NaN;
     1.11,  1.11,  NaN;
     0.23,  0.02,  0.68];

% Panel B columns: alpha (annualized %), HAC t(alpha), Sharpe ratio.
fstB = [
     2.50,  2.89,  0.54;
     6.47,  5.56,  0.94;
     5.69,  4.95,  0.85;
     2.88,  3.46,  0.71;
     4.79,  4.96,  0.85;
     5.86,  5.00,  0.86;
     6.71,  6.70,  1.05;
     8.53,  6.69,  0.99;
     2.35,  1.88,  0.46];

cfnpzB = [
    -0.79, -0.67,  0.39;
     3.68,  2.97,  0.58;
     2.65,  2.04,  0.48;
     1.01,  0.95,  0.40;
     2.28,  2.69,  0.57;
     3.85,  3.04,  0.57;
     2.39,  1.98,  0.47;
     2.83,  1.92,  0.47;
     2.35,  1.88,  0.46];

% cwMatFixed / econMatFixed pack 7 univariates + 3 combinations.
% The reply layout omits erL (row 7).
useFixed = [1:6, 8:10];
tab1.panelA = [fstA, cfnpzA, cwMatFixed(useFixed,:), reshape(cwMatAdaptive(:,:,1:2), [9,6])];
tab1.panelB = [fstB, cfnpzB, econMatFixed(useFixed,:), ...
    alphMatAdaptive(1,1:9)', tStatAlphMatAdaptive(1,1:9)', SRAdaptive(1,1:9)', ...
    alphMatAdaptive(2,1:9)', tStatAlphMatAdaptive(2,1:9)', SRAdaptive(2,1:9)'];

% Panel A: significance from the absolute CW value (one-sided test).
pval.panelA = normcdf(-abs(tab1.panelA));

% Panel B: stars only on alpha cells (cols 1, 4, 7, 10, 13), driven by
% the corresponding t(alpha) cell.
pval.panelB = ones(size(tab1.panelB));
alphaCols  = [1, 4, 7, 10, 13];
tAlphaCols = [2, 5, 8, 11, 14];
pval.panelB(:, alphaCols) = normcdf(-abs(tab1.panelB(:, tAlphaCols)));

rowLabels  = {'dp','tbl','tsp','rvar','mv','pc','comb1','comb2','comb3'};
colLabelsA = [{'Variables'}, repmat({'FS','IP','OOP'}, 1, 5)];
colLabelsB = [{'Variables'}, repmat({'$\hat{\alpha}$','$t_{\hat{\alpha}}$','Sharpe Ratio'}, 1, 5)];

LatexTableFull(tab1.panelA, colLabelsA, rowLabels, '9.2f', pval.panelA, 1)
LatexTableFull(tab1.panelB, colLabelsB, rowLabels, '9.2f', pval.panelB, 0)
end
