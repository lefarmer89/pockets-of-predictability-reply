function table4()
%TABLE4  Render Table 4: asset-pricing model OOS simulations.
%
% Compares sample CW / alpha / Sharpe statistics for the published
% sign-spec=1 forecasts against simulated distributions from five
% asset-pricing models (BY, CC, GP, W, W-no-disasters).
%
% Reads:  results/daily/OOSResults_1_2.5yE_1yDM.mat
%         results/assetPricing/<model>_asset_pricing_sims_*.mat
% Writes: stdout (LatexTableFull captured by saveTableScript).

folder = fullfile('results', 'daily');
apLabels = {'BY', 'CC', 'GP', 'W', 'W_nd'};

tab4 = NaN(21, 16);

S = load(fullfile(folder, 'OOSResults_1_2.5yE_1yDM.mat'), ...
    'statMat', 'statDiffMat', 'econMat');
refMat = [S.statMat([1 2 4], :), S.statDiffMat([1 2 4]), S.econMat([1 2 4], :)]';
tab4(:, 1) = refMat(:);

for ii = 1:5
    fileName = fullfile('results', 'assetPricing', sprintf( ...
        '%s_asset_pricing_sims_signRestriction_0_coefRestriction_0_25ybandwidth_1yDM_OOS.mat', ...
        apLabels{ii}));
    M = load(fileName, 'cwMat', 'cwDiffMat', 'econMat', 'tStatAlphaMat');

    % t(alpha) lives in tStatAlphaMat; copy into the third econMat slot
    % so the layout matches the reference (alpha, t-alpha, SR).
    econMat = M.econMat;
    econMat(:, :, 3, 1) = econMat(:, :, 2, 1);
    econMat(:, :, 2, 1) = M.tStatAlphaMat(:, :, 1);

    meanMat = [squeeze(mean(M.cwMat, 'omitnan')), mean(M.cwDiffMat, 'omitnan')', ...
        squeeze(mean(econMat(:, :, :, 1), 'omitnan'))]';
    seMat = [squeeze(std(M.cwMat, 'omitnan')), std(M.cwDiffMat, 'omitnan')', ...
        squeeze(std(econMat(:, :, :, 1), 'omitnan'))]';
    pMat = [squeeze(mean(M.cwMat > reshape(kron(refMat(1:3, :)', ...
        ones(size(M.cwMat, 1), 1)), [size(M.cwMat, 1), 3, 3]), 'omitnan')), ...
        mean(M.cwDiffMat > refMat(4, :), 'omitnan')', ...
        squeeze(mean(econMat(:, :, :, 1) > reshape(kron(refMat(5:7, :)', ...
        ones(size(econMat, 1), 1)), [size(econMat, 1), 3, 3]), 'omitnan'))]';
    tab4(:, (ii-1)*3+2:ii*3+1) = [meanMat(:), seMat(:), pMat(:)];
end

rowLabels = repmat({'$CW_{FS}$', '$CW_{IP}$', '$CW_{OOP}$', '$CW_{DIFF}$', ...
    '$\hat{\alpha}$', '$t_{\hat{\alpha}}$', 'SR'}, 1, 3);
colLabels = [{'Stats', 'Sample'}, repmat({'Avg.', 'Std. err.', 'p-val'}, 1, 5)];

LatexTableFull(tab4, colLabels, rowLabels, '9.2f', [], 0)
end
