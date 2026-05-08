function tableA4()
%TABLEA4  Render Table A.4: Fama-French SMB and HML factor portfolios.
%
% Same layout as Table 3 but with the SMB and HML target series in place
% of the market excess return.
%
% Reads:  results/famaFrench/OOSResults{SMB,HML}_1_2.5yE_1yDM.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

ffFolder = fullfile('results', 'famaFrench');

tabA4.panelA = NaN(9, 8);
tabA4.panelB = NaN(10, 8);

S = load(fullfile(ffFolder, 'OOSResultsSMB_1_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tabA4.panelA(:, 1:4) = [S.statMat(1:9, :), S.statDiffMat(1:9)];
tabA4.panelB(:, 1:3) = S.econMat([1:9 11], :);

S = load(fullfile(ffFolder, 'OOSResultsHML_1_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
tabA4.panelA(:, 5:8) = [S.statMat(1:9, :), S.statDiffMat(1:9)];
tabA4.panelB(:, 5:7) = S.econMat([1:9 11], :);

pval.panelA = normcdf(-abs(tabA4.panelA));
pval.panelB = ones(10, 8);
pval.panelB(:, [1 5]) = normcdf(-abs(tabA4.panelB(:, [2 6])));

rowLabelsA = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
rowLabelsB = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', 'pm'};
colLabelsA = {'Variables', 'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference', ...
    'Full sample', 'In-pocket', 'Out-of-pocket', 'Difference'};
colLabelsB = {'Variables', ...
    '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'Sharpe Ratio', '', ...
    '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'Sharpe Ratio', ''};

LatexTableFull(tabA4.panelA, colLabelsA, rowLabelsA, '9.2f', pval.panelA, 1)
LatexTableFull(tabA4.panelB, colLabelsB, rowLabelsB, '9.2f', pval.panelB, 1)
end
