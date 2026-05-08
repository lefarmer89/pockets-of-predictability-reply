function table3()
%TABLE3  Render Table 3: sign-restriction comparison.
%
% Panel A reports Clark-West statistics (full sample, in-pocket,
% out-of-pocket, in-out difference) across three sign-restriction specs.
% Panel B reports alpha, t(alpha), and Sharpe ratio across the same three
% specs, plus a prevailing-mean reference row.
%
% Reads:  results/daily/OOSResults_{1,2,3}_2.5yE_1yDM.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

folder = fullfile('results', 'daily');

tab3.panelA = NaN(9, 12);
tab3.panelB = NaN(10, 12);

% statMat / econMat row layout: 1-6 univariate (dp, tbl, tsp, rvar, mv,
% pc), 7 erL, 8-10 comb1/comb2/comb3, 11 lpm, 12 pm, 13 ar1. The reply
% omits erL. Panel B's 10th row is the prevailing mean (row 12).
selA = [1:6, 8:10];
selB = [1:6, 8:10, 12];

S1 = load(fullfile(folder, 'OOSResults_1_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tab3.panelA(:, 1:4) = [S1.statMat(selA, :), S1.statDiffMat(selA, :)];
tab3.panelB(:, 1:3) = S1.econMat(selB, :);

S2 = load(fullfile(folder, 'OOSResults_2_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tab3.panelA(:, 5:8) = [S2.statMat(selA, :), S2.statDiffMat(selA, :)];
tab3.panelB(:, 5:7) = S2.econMat(selB, :);

S3 = load(fullfile(folder, 'OOSResults_3_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tab3.panelA(:, 9:12) = [S3.statMat(selA, :), S3.statDiffMat(selA, :)];
tab3.panelB(:, 9:11) = S3.econMat(selB, :);

pval.panelA = normcdf(-abs(tab3.panelA));
pval.panelB = ones(10, 12);
pval.panelB(:, [1 5 9]) = normcdf(-abs(tab3.panelB(:, [2 6 10])));

rowLabelsA = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
rowLabelsB = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', 'pm'};
colLabelsA = {'Variables', 'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference', ...
    'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference', ...
    'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference'};
colLabelsB = {'Variables', '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'Sharpe Ratio', '', ...
    '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'Sharpe Ratio', '', ...
    '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'Sharpe Ratio', ''};

LatexTableFull(tab3.panelA, colLabelsA, rowLabelsA, '9.2f', pval.panelA, 1)
LatexTableFull(tab3.panelB, colLabelsB, rowLabelsB, '9.2f', pval.panelB, 0)
end
