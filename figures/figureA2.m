function figureA2()
%FIGUREA2  Render Figure A.2: hyperparameter histograms of in-pocket CW
% across the 9,720-combination grid for the four main predictors.
% Three columns per row: top 100 vs bottom 100, top 100 vs full, top 100
% vs full restricted to combos sharing the published pocket rules.
%
% Reads:  results/aggregates/topKbotK.mat
% Writes: gcf (saved to .eps by saveFigureScript).

S = load('results/aggregates/topKbotK.mat', ...
    'cwInMatTopK', 'cwInMatBotK', 'cwMat', 'cwBenchmark', 'paramCombs');

labelVec = {'dy', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
nbins = 40;

figure('Name', 'A.2')

for ii = 1:4
    d1 = [squeeze(S.cwInMatTopK(:, ii, 1)); squeeze(S.cwInMatTopK(:, ii, 2))];
    d2 = [squeeze(S.cwInMatBotK(:, ii, 1)); squeeze(S.cwInMatBotK(:, ii, 2))];
    dAll = squeeze(S.cwMat(ii, 2, :));
    dSamePocket = squeeze(S.cwMat(ii, 2, ...
        S.paramCombs(:, 6) == 1 & S.paramCombs(:, 7) == 0 & S.paramCombs(:, 8) == 21));
    refMin = @(x) 0.5 * floor(2 * min([d1(:); d2(:); x(:); ...
        S.cwBenchmark(ii, 2); squeeze(S.cwInMatTopK(1, ii, :))]));
    refMax = @(x) 0.5 * ceil( 2 * max([d1(:); d2(:); x(:); ...
        S.cwBenchmark(ii, 2); squeeze(S.cwInMatTopK(1, ii, :))]));
    drawXLines = @() arrayfun(@(f) f(), { ...
        @() xline(S.cwInMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.'), ...
        @() xline(S.cwInMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':'), ...
        @() xline(S.cwBenchmark(ii, 2),    'Color', 'k', 'LineWidth', 4)}); %#ok<NASGU>

    % --- col 1: top-100 vs bottom-100
    subplot(4, 3, 3*(ii-1)+1)
    edges = linspace(refMin(d2), refMax(d2), nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'k',   'FaceAlpha', 0.5,  'LineWidth', 1.5);
    hold on
    histogram(d2, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'red', 'FaceAlpha', 0.35, 'LineStyle', '-.');
    title([labelVec{ii}, ': In-Pocket CW Statistic'])
    xlabel('CW statistic -- in pocket'); ylabel('relative frequency')
    xline(S.cwInMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.cwInMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.cwBenchmark(ii, 2),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('top 100', 'bot 100', 'Location', 'NorthWest')
    end

    % --- col 2: top-100 vs all
    subplot(4, 3, 3*(ii-1)+2)
    edges = linspace(refMin(dAll), refMax(dAll), nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'k',    'FaceAlpha', 0.5,  'LineWidth', 1.5);
    hold on
    histogram(dAll, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'blue', 'FaceAlpha', 0.35, 'LineStyle', ':');
    title([labelVec{ii}, ': In-Pocket CW Statistic'])
    xlabel('CW statistic -- in pocket'); ylabel('relative frequency')
    xline(S.cwInMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.cwInMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.cwBenchmark(ii, 2),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('top 100', 'all', 'Location', 'NorthWest')
    end

    % --- col 3: top-100 vs all-with-same-pocket-rules
    subplot(4, 3, 3*(ii-1)+3)
    edges = linspace(refMin(dSamePocket), refMax(dSamePocket), nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'k', 'FaceAlpha', 0.5, 'LineWidth', 1.5);
    hold on
    histogram(dSamePocket, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', [0.20 0.65 0.20], 'FaceAlpha', 0.35, 'LineStyle', ':');
    title([labelVec{ii}, ': In-Pocket CW Statistic'])
    xlabel('CW statistic -- in pocket'); ylabel('relative frequency')
    xline(S.cwInMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.cwInMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.cwBenchmark(ii, 2),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('top 100', 'all w/ same pocket rules', 'Location', 'NorthWest')
    end
end
end
