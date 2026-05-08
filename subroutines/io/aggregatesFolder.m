function folder = aggregatesFolder(paths)
% AGGREGATESFOLDER  Resolve paths.results/aggregates and ensure it exists.
%
%   folder = aggregatesFolder(paths)
%
% Aggregates contains cached `.mat` files that are written by one entry-
% point and consumed by display-side renderers across multiple tables /
% figures: tab1Results, tabA3Results, topKbotK[*], oosCAlphas[*],
% signifCell, cComparisons, topKbotK_pre1989, topKbotK_post1989.

folder = fullfile(paths.results, 'aggregates');

if ~isfolder(folder)
    mkdir(folder);
end
end
