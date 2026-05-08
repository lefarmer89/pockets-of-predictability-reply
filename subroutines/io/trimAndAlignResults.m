function R = trimAndAlignResults(R, divisor, fieldsToScale)
% TRIMANDALIGNRESULTS  Post-load trimming + return-scale division for the
% empirics functions' "forecast combination" stage.
%
%   R = trimAndAlignResults(R, divisor)
%   R = trimAndAlignResults(R, divisor, fieldsToScale)
%
% Inputs:
%   R              struct of loaded forecast-result matrices. Required
%                  fields: yF2Mat (used to size the trim), and any of
%                  {fDMMat, r2Mat, r2GamMat, yF1Mat, yF1PMMat, yF2Mat,
%                  aMat, yActual, riskFree, dateVec, pocketIndMat,
%                  yAR1, rvar} that are present.
%   divisor        scalar to divide return-scaled fields by (typically 100
%                  to convert percentage returns to decimals).
%   fieldsToScale  optional cellstr of field names that should be divided
%                  by `divisor` after trimming. Defaults to
%                  {'yActual','yAR1','yF1Mat','yF1PMMat','yF2Mat','riskFree'}.
%                  Fields not present in R are skipped.
%
% Output:
%   R  same struct, with each field trimmed of `trim` leading rows
%      (where trim = max(sum(isnan(yF2Mat)))) and the scaled fields
%      divided by `divisor`. Sets pocketIndMat NaNs to 0.
%
% Replaces the trim block at e.g. dailyEmpirics.m:445-462.

if nargin < 3 || isempty(fieldsToScale)
    fieldsToScale = {'yActual','yAR1','yF1Mat','yF1PMMat','yF2Mat','riskFree'};
end

% Trim equals the longest leading-NaN column count in yF2Mat.
trim = max(sum(isnan(R.yF2Mat)));

trimFields = {'fDMMat','r2Mat','r2GamMat','yF1Mat','yF1PMMat','yF2Mat', ...
              'aMat','yActual','riskFree','dateVec','rvar','pocketIndMat','yAR1'};
for k = 1:numel(trimFields)
    f = trimFields{k};
    if isfield(R, f) && ~isempty(R.(f))
        if size(R.(f), 1) > trim
            R.(f) = R.(f)(trim+1:end, :);
        end
    end
end

% Division by divisor on the scaled return fields.
for k = 1:numel(fieldsToScale)
    f = fieldsToScale{k};
    if isfield(R, f) && ~isempty(R.(f))
        R.(f) = R.(f) ./ divisor;
    end
end

% Convention: the saved pocket indicator may have NaNs from short
% predictor histories. Treat NaN as "not in pocket".
if isfield(R, 'pocketIndMat')
    R.pocketIndMat(isnan(R.pocketIndMat)) = 0;
end
end
