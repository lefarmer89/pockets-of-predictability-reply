function tableA6()
%TABLEA6  Render Table A.6: monthly out-of-sample forecasting performance.
%
% Reports Panel A (CW full / in-pocket / out-of-pocket / difference) and
% Panel B (alpha / t(alpha) / Sharpe with prevailing-mean reference)
% across three sign-restriction specs at monthly frequency.
%
% Reads:  results/monthly/OOSResultsMonthly_{1,2,3}_2.5yE_1yDM.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

monthlyFolder = fullfile('results', 'monthly');

tabA6.panelA = NaN(9,  12);
tabA6.panelB = NaN(10, 12);

S = load(fullfile(monthlyFolder, 'OOSResultsMonthly_1_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tabA6.panelA(:, 1:4) = [S.statMat(1:9, :), S.statDiffMat(1:9, :)];
tabA6.panelB(:, 1:3) = S.econMat([1:9 11], :);

S = load(fullfile(monthlyFolder, 'OOSResultsMonthly_2_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tabA6.panelA(:, 5:8) = [S.statMat(1:9, :), S.statDiffMat(1:9, :)];
tabA6.panelB(:, 5:7) = S.econMat([1:9 11], :);

S = load(fullfile(monthlyFolder, 'OOSResultsMonthly_3_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tabA6.panelA(:, 9:12)  = [S.statMat(1:9, :), S.statDiffMat(1:9, :)];
tabA6.panelB(:, 9:11)  = S.econMat([1:9 11], :);

pval.panelA = normcdf(-abs(tabA6.panelA));
pval.panelB = ones(10, 12);
pval.panelB(:, [1 5 9]) = normcdf(-abs(tabA6.panelB(:, [2 6 10])));

rowLabelsA = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
rowLabelsB = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', 'pm'};
colLabels = {'Variables', 'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference', ...
    'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference', ...
    'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference'};

LatexTableFull(tabA6.panelA, colLabels, rowLabelsA, '9.2f', pval.panelA, 1)
LatexTableFull(tabA6.panelB, colLabels, rowLabelsB, '9.2f', pval.panelB, 0)
end
