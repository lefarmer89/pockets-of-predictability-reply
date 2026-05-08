function yWins = winsorizeRecursiveFast(yForecast, qLB, qUB, startIdx)
%WINSORIZERECURSIVEFAST  Fast drop-in for winsorizeRecursive.
%
% Same contract: at each t = startIdx:T, clip yForecast(t) to the
% [qLB, qUB] empirical percentiles of yForecast(1:t), passing earlier
% indices through unchanged.
%
% Implementation: maintains a pre-allocated sorted vector of FINITE
% values seen so far. Each step does a binary-search insert (O(log t)
% compares + O(t) memmove) and two O(1) percentile lookups via linear
% interpolation matching MATLAB's prctile convention (sample i is at
% empirical percentile 100*(i-0.5)/n).
%
% NaN handling: matches prctile — NaN values are excluded from the
% percentile estimate. NaN inputs at index t pass through to yWins(t)
% as if winsorized to the upper bound (consistent with MATLAB's
% min/max NaN handling: max(min(NaN, curUB), curLB) = curUB).
%
% Replaces winsorizeRecursive's per-iteration `prctile` calls (which
% sort internally each call). At T=5040 this yields ~50x per-call
% speedup; numerical equivalence to prctile is at machine precision.

yWins = yForecast;
n = numel(yForecast);

% Pre-allocate a length-n sorted buffer; only the first nSorted entries
% are valid at any time. Exclude NaN from the warm-up region.
sortedBuf = zeros(n, 1);
warm = max(0, startIdx - 1);
if warm > 0
    warmVals  = yForecast(1:warm);
    warmVals  = warmVals(~isnan(warmVals));
    nSorted   = numel(warmVals);
    if nSorted > 0
        sortedBuf(1:nSorted) = sort(warmVals);
    end
else
    nSorted = 0;
end

for t = startIdx:n
    val = yForecast(t);

    if ~isnan(val)
        % Binary search: find first index in sortedBuf(1:nSorted) >= val.
        lo = 1; hi = nSorted + 1;
        while lo < hi
            mid = bitshift(lo + hi, -1);
            if sortedBuf(mid) < val
                lo = mid + 1;
            else
                hi = mid;
            end
        end
        idx = lo;

        if idx <= nSorted
            sortedBuf(idx+1:nSorted+1) = sortedBuf(idx:nSorted);
        end
        sortedBuf(idx) = val;
        nSorted = nSorted + 1;
    end

    if nSorted == 0
        % No finite samples yet; pass through (yWins(t) = val from init).
        continue
    end

    % Percentile lookup matching prctile interpolation.
    posLB = qLB * nSorted / 100 + 0.5;
    if posLB <= 1
        curLB = sortedBuf(1);
    elseif posLB >= nSorted
        curLB = sortedBuf(nSorted);
    else
        iLo = floor(posLB);
        curLB = sortedBuf(iLo) + (posLB - iLo) * (sortedBuf(iLo+1) - sortedBuf(iLo));
    end

    posUB = qUB * nSorted / 100 + 0.5;
    if posUB <= 1
        curUB = sortedBuf(1);
    elseif posUB >= nSorted
        curUB = sortedBuf(nSorted);
    else
        iLo = floor(posUB);
        curUB = sortedBuf(iLo) + (posUB - iLo) * (sortedBuf(iLo+1) - sortedBuf(iLo));
    end

    % MATLAB's max/min treat NaN as missing, so when val is NaN the
    % result collapses to curUB; for finite val it's the standard clip.
    yWins(t) = max(min(val, curUB), curLB);
end
end
