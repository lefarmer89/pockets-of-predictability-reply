function table2()
%TABLE2  Render Table 2: pocket statistics (daily).
%
% Columns 1-4 hold all-pockets statistics (count, in-sample fraction,
% duration min/mean/max, integral R^2 min/mean/max). Columns 5-8 restrict
% to pockets flagged as significant by the dailyBootstrap-driven
% signifCell.
%
% Reads:  results/daily/forecastResults_1_2.5yE_1yDM_1S_pm.mat
%         results/aggregates/signifCell.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

folder = fullfile('results', 'daily');
load(fullfile(folder, 'forecastResults_1_2.5yE_1yDM_1S_pm.mat'), ...
    'durationMat', 'integralR2Mat', 'pocketIndMat');

% Add the 21-day minimum window back into reported duration.
durationMat   = durationMat + 21;
durationMat   = durationMat(:, 1:4);
integralR2Mat = integralR2Mat(:, 1:4);
pocketIndMat  = pocketIndMat(:, 1:4);

tab2 = NaN(8);
tab2(:, 1:4) = [
    sum(~isnan(durationMat));
    mean(pocketIndMat, 'omitnan');
    min(durationMat, [], 'omitnan');
    mean(durationMat, 'omitnan');
    max(durationMat, [], 'omitnan');
    min(integralR2Mat, [], 'omitnan');
    mean(integralR2Mat, 'omitnan');
    max(integralR2Mat, [], 'omitnan')];

load(fullfile('results', 'aggregates', 'signifCell.mat'), 'signifCell');
% signifCell uses the internal duration that excludes the 21-day minimum,
% so subtract here to align indexing.
durationMat = durationMat - 21;

for ii = 1:4
    keepInd = logical(signifCell{ii});
    keepInd2 = zeros(size(pocketIndMat, 1), 1);
    if pocketIndMat(1, ii) == 1
        curPocket = 1;
        inPocket  = true;
    else
        curPocket = 0;
        inPocket  = false;
    end

    for t = 2:size(pocketIndMat, 1)
        if (pocketIndMat(t, ii) - pocketIndMat(t-1, ii)) == 1
            curPocket = curPocket + 1;
            inPocket  = true;
        end
        if (pocketIndMat(t, ii) - pocketIndMat(t-1, ii)) == -1
            inPocket = false;
        end
        if inPocket && (keepInd(curPocket) == 1)
            keepInd2(t) = 1;
        end
    end
    keepInd2(isnan(pocketIndMat(:, ii))) = NaN;
    tab2(:, ii+4) = [
        sum(~isnan(durationMat(keepInd, ii)));
        mean(keepInd2, 'omitnan');
        min(durationMat(keepInd, ii), [], 'omitnan');
        mean(durationMat(keepInd, ii), 'omitnan');
        max(durationMat(keepInd, ii), [], 'omitnan');
        min(integralR2Mat(keepInd, ii), [], 'omitnan');
        mean(integralR2Mat(keepInd, ii), 'omitnan');
        max(integralR2Mat(keepInd, ii), [], 'omitnan')];
end

rowLabels = {'Num pockets', 'Fraction of sample', ...
    '\hspace{1em} Min', '\hspace{1em} Mean', '\hspace{1em} Max', ...
    '\hspace{1em} Min', '\hspace{1em} Mean', '\hspace{1em} Max'};
colLabels = {'Statistics', 'dp', 'tbl', 'tsp', 'rvar', 'dp', 'tbl', 'tsp', 'rvar'};

LatexTableFull(tab2, colLabels, rowLabels, '9.2f', [], 0)
end
