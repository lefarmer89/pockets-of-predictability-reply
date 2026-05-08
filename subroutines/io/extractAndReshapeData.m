function S = extractAndReshapeData(rawData, TMax, B)
% EXTRACTANDRESHAPEDATA  Reshape the sticky-expectations 5-column raw
% simulation matrix into separate 3D arrays for excess return, true
% expected return, and the predictor stack.
%
%   S = extractAndReshapeData(rawData, TMax, B)
%
% Inputs:
%   rawData  Raw matrix loaded from the calibration .mat. Layout follows
%            the package convention: TMax rows per replication, B
%            replications stacked vertically; columns 1..5 hold the
%            simulated state variables in the canonical order.
%   TMax     observations per replication
%   B        number of replications
%
% Output struct S with fields:
%   exretMat       TMax x B excess return
%   trueERMat      TMax x B true expected return
%   predictorMat   TMax x B x 3 stacked predictors (dp, rf, rvar)

S.exretMat     = reshape(rawData(:, 1), TMax, B);
S.trueERMat    = reshape(rawData(:, 2), TMax, B);

predictorRaw   = rawData(:, 3:5);
predictorMat   = NaN(TMax, B, 3);
for k = 1:3
    predictorMat(:, :, k) = reshape(predictorRaw(:, k), TMax, B);
end
S.predictorMat = predictorMat;
end
