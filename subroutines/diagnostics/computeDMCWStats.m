function row = computeDMCWStats(diffSeries, pocketIdx, kind)
% COMPUTEDMCWSTATS  Per-row HAC t-stat triple for one forecast differential.
%
%   row = computeDMCWStats(diffSeries, pocketIdx, kind)
%
% Inputs:
%   diffSeries  N x 1 forecast differential (DM or CW).
%   pocketIdx   N x 1 logical / 0-1 pocket indicator. Pass [] or all-true
%               to skip pocket segmentation (full-sample only).
%   kind        'dm' returns 1 x 3 [full, in-pocket, out-of-pocket]
%               'cw' returns the same plus the Welch-difference t-stat
%                    appended as the 4th element.
%
% Output:
%   row  1 x 3 (kind = 'dm') or 1 x 4 (kind = 'cw') of HAC t-stats.
%        Empty subsamples emit NaN.
%
% Used by the bootstrap and hyperparameter functions where the caller
% loops over replications or combos rather than predictor specs.

if isempty(pocketIdx)
    pocketIdx = true(size(diffSeries));
end
pIn  = logical(pocketIdx);
pOut = ~pIn;
nIn  = sum(pIn);
nOut = sum(pOut);

% Full sample
[tFull, ~, ~] = hacTwithSE(diffSeries);

% In-pocket
if nIn > 0
    [tIn, seIn, mIn] = hacTwithSE(diffSeries(pIn));
else
    tIn = NaN; seIn = NaN; mIn = NaN;
end

% Out-of-pocket
if nOut > 0
    [tOut, seOut, mOut] = hacTwithSE(diffSeries(pOut));
else
    tOut = NaN; seOut = NaN; mOut = NaN;
end

switch lower(kind)
    case 'dm'
        row = [tFull, tIn, tOut];
    case 'cw'
        if nIn > 1 && nOut > 1 && ~isnan(seIn) && ~isnan(seOut)
            seDiff = sqrt(((nIn-1)*seIn^2 + (nOut-1)*seOut^2) / (nIn+nOut-2));
            tDiff  = (mIn - mOut) / (seDiff * sqrt(1/nIn + 1/nOut));
        else
            tDiff = NaN;
        end
        row = [tFull, tIn, tOut, tDiff];
    otherwise
        error('computeDMCWStats:badKind', 'kind must be ''dm'' or ''cw''');
end
end


function [t, se, m] = hacTwithSE(x)
if isempty(x)
    t = NaN; se = NaN; m = NaN;
    return
end
s  = regstats2Fast(x, ones(numel(x), 1), 'onlydata', 'hac');
% regstats2's hac.t/hac.se are 2-element vectors only when p>=2; for the
% p=1 'onlydata' design used here both Fast and reference return scalars,
% so no indexing change needed.
t  = s.hac.t;
se = s.hac.se * sqrt(numel(x) - 1);
m  = mean(x, 'omitmissing');
end
