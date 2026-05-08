function figureA1()
%FIGUREA1  Render Figure A.1: hyperparameter histograms of t(alpha)
% across the 9,720-combination grid for the four main predictors.
% Left column: top 100 vs bottom 100 (both selectors stacked). Right
% column: top 100 vs full distribution.
%
% Reads:  results/aggregates/topKbotK.mat
% Writes: gcf (saved to .eps by saveFigureScript).

S = load('results/aggregates/topKbotK.mat', ...
    'tStatAlphMatTopK', 'tStatAlphMatBotK', 'tStatAlphMatEnd', ...
    'tAlphAllSpecs', 'tStatAlphaBenchmark');

labelVec = {'dy', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};
nbins = 40;

figure('Name', 'A.1')

for ii = 1:4
    % --- Left: top 100 vs bottom 100
    subplot(4, 2, 2*(ii-1)+1)
    d1 = [squeeze(S.tStatAlphMatTopK(:, ii, 1)); squeeze(S.tStatAlphMatTopK(:, ii, 2))];
    d2 = [squeeze(S.tStatAlphMatBotK(:, ii, 1)); squeeze(S.tStatAlphMatBotK(:, ii, 2))];
    d3 = S.tStatAlphMatEnd(ii, :);
    edges = linspace( ...
        0.5 * floor(2 * min([d1(:); d2(:); d3(:); S.tStatAlphaBenchmark(ii); squeeze(S.tStatAlphMatTopK(1, ii, :))])), ...
        0.5 * ceil( 2 * max([d1(:); d2(:); d3(:); S.tStatAlphaBenchmark(ii); squeeze(S.tStatAlphMatTopK(1, ii, :))])), ...
        nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'k',   'FaceAlpha', 0.5,  'LineWidth', 1.5);
    hold on
    histogram(d2, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'red', 'FaceAlpha', 0.35, 'LineStyle', '-.');
    title([labelVec{ii}, ': t-stat of \alpha'])
    xlabel('t stat \alpha')
    ylabel('relative frequency')
    xline(S.tStatAlphMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.tStatAlphMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.tStatAlphaBenchmark(ii),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('top 100', 'bot 100', 'Location', 'NorthWest')
    end

    % --- Right: top 100 vs full distribution
    subplot(4, 2, 2*(ii-1)+2)
    d1 = [squeeze(S.tStatAlphMatTopK(:, ii, 1)); squeeze(S.tStatAlphMatTopK(:, ii, 2))];
    d2 = [squeeze(S.tStatAlphMatBotK(:, ii, 1)); squeeze(S.tStatAlphMatBotK(:, ii, 2))]; %#ok<NASGU>
    d3 = S.tAlphAllSpecs(:, ii);
    edges = linspace( ...
        0.5 * floor(2 * min([d1(:); d3(:); S.tStatAlphaBenchmark(ii); squeeze(S.tStatAlphMatTopK(1, ii, :))])), ...
        0.5 * ceil( 2 * max([d1(:); d3(:); S.tStatAlphaBenchmark(ii); squeeze(S.tStatAlphMatTopK(1, ii, :))])), ...
        nbins+1);
    histogram(d1, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'k',    'FaceAlpha', 0.5,  'LineWidth', 1.5);
    hold on
    histogram(d3, 'Normalization', 'probability', 'BinEdges', edges, ...
        'FaceColor', 'blue', 'FaceAlpha', 0.35, 'LineStyle', ':');
    title([labelVec{ii}, ': t-stat of \alpha'])
    xlabel('t stat \alpha')
    ylabel('relative frequency')
    xline(S.tStatAlphMatTopK(1, ii, 1), 'blue', 'LineWidth', 3, 'LineStyle', '-.')
    xline(S.tStatAlphMatTopK(1, ii, 2), 'red',  'LineWidth', 3, 'LineStyle', ':')
    xline(S.tStatAlphaBenchmark(ii),    'Color', 'k', 'LineWidth', 4)
    xlim('tight')
    if ii == 1
        legend('top 100', 'all', 'Location', 'NorthWest')
    end
end
end
