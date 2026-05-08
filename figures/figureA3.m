function figureA3()
%FIGUREA3  Render Figure A.3: scatter plot of t(alpha) under true OOS
% c-scaling vs ex-post c-scaling for the four main predictors. Diagonal
% indicates exact agreement.
%
% Reads:  results/aggregates/topKbotK.mat
% Writes: gcf (saved to .eps by saveFigureScript).

S = load('results/aggregates/topKbotK.mat', ...
    'tAlphAllSpecs', 'tStatAlphMatEnd');

labelVec = {'dy', 'tbl', 'tsp', 'rvar', 'mv', 'pc', 'comb1', 'comb2', 'comb3'};

figure('Name', 'A.3')
for ii = 1:4
    subplot(2, 2, ii)
    scatter(S.tAlphAllSpecs(:, ii), squeeze(S.tStatAlphMatEnd(ii, :)))
    hold on
    rangetmp = linspace( ...
        min([S.tAlphAllSpecs(:, ii); S.tStatAlphMatEnd(ii, :)']), ...
        max([S.tAlphAllSpecs(:, ii); S.tStatAlphMatEnd(ii, :)']), 50)';
    plot(rangetmp, rangetmp, 'red')
    title([labelVec{ii}, ': t Stats on alpha'])
    xlabel('true OOS')
    ylabel('ex post c')
end
end
