function tableA1()
%TABLEA1  Render Table A.1: G-prior parameter robustness (g in {1, 2, 3}).
%
% Stacks Panel A (Clark-West full / in-pocket / out-of-pocket) and Panel B
% (alpha / t(alpha) / Sharpe) across three g values for the published
% sign-spec=1 / 2.5y kernel / 1y SED baseline.
%
% Reads:  results/gPriorRobustness/{g1,g3}/OOSResults_2_2.5yE_1yDM.mat
%         results/daily/OOSResults_2_2.5yE_1yDM.mat   (g=2 baseline)
% Writes: stdout (LatexTableFull captured by saveTableScript).

folders = {fullfile('results', 'gPriorRobustness', 'g1'), ...
           fullfile('results', 'daily'), ...
           fullfile('results', 'gPriorRobustness', 'g3')};

tabA1 = NaN(18, 9);

% Daily statMat / econMat row layout: 1-6 univariate, 7 erL, 8-10 combs.
% The reply table omits erL.
selA = [1:6, 8:10];

for jj = 1:3
    fileName = fullfile(folders{jj}, 'OOSResults_2_2.5yE_1yDM.mat');
    S = load(fileName, 'statMat', 'statDiffMat', 'econMat');
    tabA1(1:9,    jj)   = S.statMat(selA, 2);
    tabA1(1:9,    jj+3) = S.statMat(selA, 3);
    tabA1(1:9,    jj+6) = S.statDiffMat(selA);
    tabA1(10:18,  jj)   = S.econMat(selA, 1);
    tabA1(10:18,  jj+3) = S.econMat(selA, 2);
    tabA1(10:18,  jj+6) = S.econMat(selA, 3);
end

bottomHalf = NaN(9, 9);
bottomHalf(:, 1) = normcdf(-abs(tabA1(10:18, 4)));
bottomHalf(:, 2) = normcdf(-abs(tabA1(10:18, 5)));
bottomHalf(:, 3) = normcdf(-abs(tabA1(10:18, 6)));

pval = [normcdf(-abs(tabA1(1:9, :))); bottomHalf];

rowLabels = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', ...
             'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
colLabels = [{'Variables'}, repmat({'$g=1$', '$g=2$', '$g=3$'}, 1, 3)];

LatexTableFull(tabA1, colLabels, rowLabels, '9.2f', pval, 1)
end
