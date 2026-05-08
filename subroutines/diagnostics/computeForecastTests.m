function [dmMat, cwMat, cwDiffMat] = computeForecastTests(fDMCell, fCWCell, pocketIndCell)
% COMPUTEFORECASTTESTS  HAC t-stats for forecast differentials,
% computed for the full sample, in-pocket, and out-of-pocket subsamples.
%
%   [dmMat, cwMat, cwDiffMat] = computeForecastTests(fDMCell, fCWCell, pocketIndCell)
%
% Inputs (each {n,1} cell array, all with same n):
%   fDMCell        Diebold-Mariano differential per spec
%   fCWCell        Clark-West differential per spec
%   pocketIndCell  pocket indicator (logical or 0/1) per spec
%
% Outputs:
%   dmMat      n x 3 HAC t-stat: [full, in-pocket, out-of-pocket]
%   cwMat      n x 3 HAC t-stat for CW differentials, same columns
%   cwDiffMat  n x 1 Welch-style t-stat comparing in-pocket vs
%              out-of-pocket CW means. NaN when either subset is empty.
%
% Replaces the ~300-line block of regstats2 calls in dailyEmpirics,
% monthlyEmpirics, and famaFrenchEmpirics. Empty subsets leave the
% corresponding column NaN.

n = numel(fDMCell);
dmMat     = NaN(n, 3);
cwMat     = NaN(n, 3);
cwDiffMat = NaN(n, 1);

for i = 1:n
    fDM = fDMCell{i};
    fCW = fCWCell{i};
    pIn = logical(pocketIndCell{i});
    pOut = ~pIn;

    % Full sample
    dmMat(i, 1) = hacT(fDM);
    cwMat(i, 1) = hacT(fCW);

    % In-pocket
    nIn = sum(pIn);
    if nIn > 0
        dmMat(i, 2)             = hacT(fDM(pIn));
        [cwMat(i, 2), seIn, mIn] = hacTwithSE(fCW(pIn));
    else
        seIn = NaN; mIn = NaN;
    end

    % Out-of-pocket
    nOut = sum(pOut);
    if nOut > 0
        dmMat(i, 3)               = hacT(fDM(pOut));
        [cwMat(i, 3), seOut, mOut] = hacTwithSE(fCW(pOut));
    else
        seOut = NaN; mOut = NaN;
    end

    % Welch-style difference of in-pocket vs out-of-pocket CW means.
    if nIn > 1 && nOut > 1 && ~isnan(seIn) && ~isnan(seOut)
        seDiff = sqrt(((nIn-1)*seIn^2 + (nOut-1)*seOut^2) / (nIn+nOut-2));
        cwDiffMat(i) = (mIn - mOut) / (seDiff * sqrt(1/nIn + 1/nOut));
    end
end
end


function t = hacT(x)
% HAC t-stat of x against a constant.
if isempty(x)
    t = NaN;
    return
end
s = regstats2Fast(x, ones(numel(x), 1), 'onlydata', 'hac');
t = s.hac.t;
end


function [t, se, m] = hacTwithSE(x)
% HAC t-stat plus se and sample mean (for the Welch difference).
if isempty(x)
    t = NaN; se = NaN; m = NaN;
    return
end
s  = regstats2Fast(x, ones(numel(x), 1), 'onlydata', 'hac');
t  = s.hac.t;
se = s.hac.se * sqrt(numel(x) - 1);
m  = mean(x, 'omitmissing');
end
