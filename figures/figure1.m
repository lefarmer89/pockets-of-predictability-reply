function figure1()
%FIGURE1  Render Figure 1: distribution of t(alpha) and in-pocket CW
% across the top-100 adaptive-shrinkage models for the four main
% predictors (dp, tbl, tsp, rvar). Two columns per row: t(alpha) on the
% left, in-pocket CW on the right. Blue = min RMSE selector, red = max
% alpha selector, black vertical = published-baseline benchmark.
%
% Reads:  results/aggregates/topKbotK.mat
% Writes: gcf (saved to .eps by saveFigureScript).

S = load('results/aggregates/topKbotK.mat', ...
    'tStatAlphMatTopK', 'cwInMatTopK', 'tStatAlphaBenchmark', 'cwBenchmark');

labelVec = {'dy', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};

figure(1)
nbins = 15;

for ii = 1:4
    % --- t(alpha) panel
    subplot(4, 2, 2*(ii-1)+1)
    d1 = squeeze(S.tStatAlphMatTopK(:, ii, 1));
    d2 = squeeze(S.tStatAlphMatTopK(:, ii, 2));
    edges = linspace( ...
        0.5 * floor(2 * min([d1(:); d2(:); S.tStatAlphaBenchmark(ii)])), ...
        0.5 * ceil( 2 * max([d1(:); d2(:); S.tStatAlphaBenchmark(ii)])), ...
        nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'blue', 'FaceAlpha', 0.35);
    hold on
    histogram(d2, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'red',  'FaceAlpha', 0.35);
    title([labelVec{ii}, ': t-stat of \alpha'])
    xlabel('t stat \alpha')
    ylabel('relative frequency')
    xline(S.tStatAlphMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.tStatAlphMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.tStatAlphaBenchmark(ii),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('min RMSE', 'max \alpha', 'Location', 'NorthWest')
    end

    % --- in-pocket CW panel
    subplot(4, 2, 2*(ii-1)+2)
    d1 = squeeze(S.cwInMatTopK(:, ii, 1));
    d2 = squeeze(S.cwInMatTopK(:, ii, 2));
    edges = linspace( ...
        0.5 * floor(2 * min([d1(:); d2(:); S.cwBenchmark(ii, 2)])), ...
        0.5 * ceil( 2 * max([d1(:); d2(:); S.cwBenchmark(ii, 2)])), ...
        nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'blue', 'FaceAlpha', 0.35);
    hold on
    histogram(d2, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'red',  'FaceAlpha', 0.35);
    title([labelVec{ii}, ': In-Pocket CW statistic'])
    xlabel('CW statistic -- in pocket')
    ylabel('relative frequency')
    xline(S.cwInMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.cwInMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.cwBenchmark(ii, 2),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
end
end
