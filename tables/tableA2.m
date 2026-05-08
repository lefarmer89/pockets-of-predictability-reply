function tableA2()
%TABLEA2  Render Table A.2: window-length robustness.
%
% Reports Panel A (in-pocket / out-of-pocket CW) and Panel B (CW
% difference and economic stats) across five (coefWindow, sedWindow)
% combinations: (2.5y, 1y), (2y, 1y), (3y, 1y), (2.5y, 1.5y), (2.5y, 6m).
%
% Reads:  results/daily/OOSResults_2_<window>.mat  (5 specs)
% Writes: stdout (LatexTableFull captured by saveTableScript).

folder = fullfile('results', 'daily');
rLabels  = {'_2'};
eLabels  = {'_2.5yE', '_2yE', '_3yE', '_2.5yE', '_2.5yE'};
dmLabels = {'_1yDM',  '_1yDM', '_1yDM', '_6mDM', '_1.5yDM'};

tabA2 = NaN(18, 10);
selA  = [1:6, 8:10];   % drop erL row

tAlphaMat = NaN(9, 5);
for jj = 1:5
    fileName = fullfile(folder, sprintf('OOSResults%s%s%s.mat', ...
        rLabels{1}, eLabels{jj}, dmLabels{jj}));
    S = load(fileName, 'statMat', 'statDiffMat', 'econMat');
    tabA2(1:9,    jj)   = S.statMat(selA, 2);
    tabA2(1:9,    jj+5) = S.statMat(selA, 3);
    tabA2(10:18,  jj)   = S.statDiffMat(selA);
    tabA2(10:18,  jj+5) = S.econMat(selA, 1);
    tAlphaMat(:,  jj)   = S.econMat(selA, 2);
end

pval = normcdf(-abs(tabA2));
pval(10:18, 6:end) = normcdf(-abs(tAlphaMat));

rowLabels = {'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3', ...
             'dp', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
colLabels = [{'Variables'}, repmat({'2.5yCoef, 1ySED', ...
    '2yCoef, 1ySED', '3yCoef, 1ySED', '2.5yCoef, 1.5ySED', '2.5yCoef, 6mSED'}, 1, 2)];

LatexTableFull(tabA2, colLabels, rowLabels, '9.2f', pval, 1)
end
