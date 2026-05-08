function tableA5()
%TABLEA5  Render Table A.5: monthly pocket statistics.
%
% Reports per-variable pocket count, in-sample fraction, duration
% (min / mean / max in trading days), and integral R^2 (min / mean / max).
%
% Reads:  results/monthly/forecastResultsMonthly_1_2.5yE_1yDM_1S_pm.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

monthlyFolder = fullfile('results', 'monthly');
S = load(fullfile(monthlyFolder, 'forecastResultsMonthly_1_2.5yE_1yDM_1S_pm.mat'), ...
    'durationMat', 'integralR2Mat', 'pocketIndMat');

% Convert to days (21 per month) and add the 1-month minimum back in.
durationMat   = S.durationMat(:, 1:4) * 21 + 21;
integralR2Mat = S.integralR2Mat(:, 1:4);
pocketIndMat  = S.pocketIndMat(:, 1:4);

tabA5 = [
    sum(~isnan(durationMat));
    mean(pocketIndMat, 'omitnan');
    min(durationMat, [], 'omitnan');
    mean(durationMat, 'omitnan');
    max(durationMat, [], 'omitnan');
    min(integralR2Mat, [], 'omitnan');
    mean(integralR2Mat, 'omitnan');
    max(integralR2Mat, [], 'omitnan')];

rowLabels = {'Num pockets', 'Fraction of sample', ...
    '\hspace{1em} Min', '\hspace{1em} Mean', '\hspace{1em} Max', ...
    '\hspace{1em} Min', '\hspace{1em} Mean', '\hspace{1em} Max'};
colLabels = {'Statistics', 'dp', 'tbl', 'tsp', 'rvar'};

LatexTableFull(tabA5, colLabels, rowLabels, '9.2f', [], 0)
end
