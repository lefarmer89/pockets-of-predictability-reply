function table5()
%TABLE5  Render Table 5: sticky-expectations model simulations.
%
% Compares the sample CW / alpha / Sharpe statistics for the dp predictor
% against simulated distributions from three models: SE (sticky
% expectations), RE (rational expectations), and RE_Recalibrated.
%
% Reads:  results/assetPricing/{SE,RE,RE_Recalibrated}_asset_pricing_sims_*.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

dataMoments_2Point5Y = [0.32, 1.25, -1.03, 2.00, 2.74, 1.94, 0.48, ...
    1.34, 2.29, -1.05, 3.14, 5.48, 4.01, 0.71, ...
    0.57, 1.39, -0.83, 2.07, 3.41, 2.43, 0.52]';

tab5 = NaN(21, 10);
tab5(:, 1) = dataMoments_2Point5Y;

filePrefixes = {'SE', 'RE', 'RE_Recalibrated'};

for ii = 1:numel(filePrefixes)
    fileName = fullfile('results', 'assetPricing', sprintf( ...
        '%s_asset_pricing_sims_signRestriction_0_coefRestriction_0_25ybandwidth_1yDM_OOS.mat', ...
        filePrefixes{ii}));
    M = load(fileName, 'cwMat', 'cwDiffMat', 'econMat', 'tStatAlphaMat');

    modelMoments = NaN(21, 3);
    for jj = 1:3
        tempCW     = squeeze(M.cwMat(:, jj, :));
        tempCWDiff = M.cwDiffMat(:, jj);
        tempEcon   = squeeze(M.econMat(:, jj, :, 1));
        tempTStat  = squeeze(M.tStatAlphaMat(:, jj, 1));
        tempEcon(:, 3) = tempEcon(:, 2);
        tempEcon(:, 2) = tempTStat;

        modelMoments((jj-1)*7+1:jj*7, :) = [
            [mean(tempCW), mean(tempCWDiff), mean(tempEcon)]', ...
            [std(tempCW),  std(tempCWDiff),  std(tempEcon)]', ...
            mean([tempCW, tempCWDiff, tempEcon] > dataMoments_2Point5Y((jj-1)*7+1:jj*7)')'];
    end

    tab5(:, (ii-1)*3+2:ii*3+1) = modelMoments;
end

rowLabels = repmat({'$CW_{fs}$', '$CW_{ip}$', '$CW_{oop}$', '$CW_{diff}$', ...
    '$\alpha$', '$t_{\alpha}$', 'SR'}, 1, 3);
colLabels = [{'Data'}, repmat({'Avg', 'Std. err', 'p-val'}, 1, 3)];

LatexTableFull(tab5, colLabels, rowLabels, '9.2f', [], 1)
end
