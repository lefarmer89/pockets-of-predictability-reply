function [pocketIndices, pocketPeriods, pocketLengths] = ...
    detectPocketsFromSED(fDMHat, totalLen, minDuration)
% DETECTPOCKETSFROMSED  From the smoothed SED estimate, identify runs
% of positive periods longer than minDuration and return their indicator
% vector, start/end indices, and durations. Replaces the duplicated
% pocket-detection block that appeared inside every empirics function.
%
%   [pocketIndices, pocketPeriods, pocketLengths] = ...
%       detectPocketsFromSED(fDMHat, totalLen, minDuration)
%
% Inputs:
%   fDMHat       smoothed SED estimate from lps1_v2 (the 6th output).
%                Caller is expected to have already run the kernel
%                regression that produced this.
%   totalLen     length of the underlying time series the indicator
%                should match (typically numel(r2Local)).
%   minDuration  minimum number of consecutive positive periods required
%                to retain a pocket. Pockets are dated at minDuration
%                periods after the first positive period of the run.
%
% Outputs:
%   pocketIndices   totalLen x 1 logical indicator vector.
%   pocketPeriods   N x 2 [startIdx, endIdx] of each retained pocket.
%   pocketLengths   N x 1 duration in periods of each retained pocket.

prePad = totalLen - numel(fDMHat) + 1;
hits   = find([zeros(prePad, 1); fDMHat(1:end-1) > 0]);

if isempty(hits)
    pocketIndices = false(totalLen, 1);
    pocketPeriods = zeros(0, 2);
    pocketLengths = zeros(0, 1);
    return
end

% Group consecutive hits into [startIdx, endIdx] runs.
breaks = [true; diff(hits) ~= 1];
runStarts = hits(breaks);
runEnds   = hits([find(breaks(2:end)); numel(hits)]);
periods = [runStarts, runEnds];

% Drop runs shorter than minDuration; advance the start by minDuration so
% the pocket is dated at the point the SED has been positive for
% minDuration periods.
keep = (periods(:,2) - periods(:,1) + 1) > minDuration;
periods = periods(keep, :);
periods(:,1) = periods(:,1) + minDuration;

pocketLengths = periods(:,2) - periods(:,1) + 1;
pocketPeriods = periods;

pocketIndices = false(totalLen, 1);
for k = 1:size(periods, 1)
    pocketIndices(periods(k,1):periods(k,2)) = true;
end
end
