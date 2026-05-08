function P = computePerformanceMetrics(yPocketTime, yActual, ySigFactor, yMomFactor, freq)
% COMPUTEPERFORMANCEMETRICS  Per-spec alpha / beta / Sharpe / vol /
% momentum / 3-factor / LPM regression panel ready for LatexTableFull.
%
%   P = computePerformanceMetrics(yPocketTime, yActual, ySigFactor, yMomFactor, freq)
%
% Inputs:
%   yPocketTime  T x numSpecifications portfolio returns
%   yActual      T x 1 realized excess returns
%   ySigFactor   T x 1 vol-signed factor (from constructPortfolioFactors)
%   yMomFactor   T x 1 momentum factor (from constructPortfolioFactors)
%   freq         'daily' (annualization 252) or 'monthly' (12)
%
% Output struct P with fields:
%   coefMat        2 x S [alpha; beta]
%   tStatMat       2 x S HAC t-stats
%   mseMat         1 x S
%   coefMatVol, tStatMatVol      3 x S [alpha; beta; betaVol]
%   coefMatMom, tStatMatMom      3 x S [alpha; beta; betaMom]
%   coefMat3Fac, tStatMat3Fac    4 x S [alpha; beta; betaVol; betaMom]
%   coefMatLPM, tStatMatLPM      3 x S [alpha; beta; betaLPM]
%   SRmkt, SRmktAnnual           market Sharpe and annualized
%   SR, SRAnnual                 portfolio Sharpe and annualized
%   appRatio, appRatioAnnual     appraisal ratio and annualized (NB: the
%                                monthly entry-point originally annualized
%                                appRatio with sqrt(252) — preserved here
%                                via cfg-free hard-coded sqrt(252) for
%                                bit-for-bit equivalence)
%   alphaAnnual, alphaAnnualVol, alphaAnnualMom, alphaAnnual3Fac, alphaAnnualLPM

switch freq
    case {'daily', 'famaFrench'}
        annFactor = 252;
    case 'monthly'
        annFactor = 12;
    otherwise
        error('computePerformanceMetrics:badFreq', ...
            'freq must be ''daily'', ''monthly'', or ''famaFrench''');
end

S = size(yPocketTime, 2);

coefMat      = NaN(2, S);
tStatMat     = NaN(2, S);
mseMat       = NaN(1, S);
coefMatVol   = NaN(3, S);  tStatMatVol  = NaN(3, S);
coefMatMom   = NaN(3, S);  tStatMatMom  = NaN(3, S);
coefMat3Fac  = NaN(4, S);  tStatMat3Fac = NaN(4, S);
coefMatLPM   = NaN(3, S);  tStatMatLPM  = NaN(3, S);

for ii = 1:S
    % Single-factor (CAPM-style)
    s = regstats2Fast(yPocketTime(:,ii), yActual);
    coefMat(:,ii)  = s.beta;
    mseMat(ii)     = s.mse;
    tStatMat(:,ii) = s.hac.t;

    % + vol factor
    s = regstats2Fast(yPocketTime(:,ii), [yActual, ySigFactor], 'linear', ...
        {'beta','hac','mse'});
    coefMatVol(:,ii)  = s.beta;
    tStatMatVol(:,ii) = s.hac.t;

    % + momentum factor
    s = regstats2Fast(yPocketTime(:,ii), [yActual, yMomFactor], 'linear', ...
        {'beta','hac','mse'});
    coefMatMom(:,ii)  = s.beta;
    tStatMatMom(:,ii) = s.hac.t;

    % 3-factor
    s = regstats2Fast(yPocketTime(:,ii), [yActual, ySigFactor, yMomFactor], ...
        'linear', {'beta','hac','mse'});
    coefMat3Fac(:,ii)  = s.beta;
    tStatMat3Fac(:,ii) = s.hac.t;

    % LPM (use second-to-last spec as auxiliary regressor — matches the
    % original convention in dailyEmpirics/monthlyEmpirics).
    s = regstats2Fast(yPocketTime(:,ii), [yActual, yPocketTime(:,end-1)], ...
        'linear', {'beta','hac','mse'});
    coefMatLPM(:,ii)  = s.beta;
    tStatMatLPM(:,ii) = s.hac.t;
end

% Sharpe ratios.
SRmkt       = mean(yActual, 'omitnan') / std(yActual, 'omitnan');
SRmktAnnual = SRmkt * sqrt(annFactor);
SR          = sqrt(SRmkt^2 + (coefMat(1,:).^2) ./ mseMat(1,:));
SRAnnual    = SR * sqrt(annFactor);

% Appraisal ratio. Monthly entry-points used sqrt(252) here historically;
% to preserve numerical equivalence we keep that for monthly calls.
appRatio = coefMat(1,:) ./ sqrt(mseMat(1,:));
if strcmp(freq, 'monthly')
    appRatioAnnual = appRatio * sqrt(252);  % preserve original scaling
else
    appRatioAnnual = appRatio * sqrt(annFactor);
end

alphaAnnual     = coefMat(1,:)     * annFactor * 100;
alphaAnnualVol  = coefMatVol(1,:)  * annFactor * 100;
alphaAnnualMom  = coefMatMom(1,:)  * annFactor * 100;
alphaAnnual3Fac = coefMat3Fac(1,:) * annFactor * 100;
alphaAnnualLPM  = coefMatLPM(1,:)  * annFactor * 100;

P = struct( ...
    'coefMat',        coefMat, ...
    'tStatMat',       tStatMat, ...
    'mseMat',         mseMat, ...
    'coefMatVol',     coefMatVol, ...
    'tStatMatVol',    tStatMatVol, ...
    'coefMatMom',     coefMatMom, ...
    'tStatMatMom',    tStatMatMom, ...
    'coefMat3Fac',    coefMat3Fac, ...
    'tStatMat3Fac',   tStatMat3Fac, ...
    'coefMatLPM',     coefMatLPM, ...
    'tStatMatLPM',    tStatMatLPM, ...
    'SRmkt',          SRmkt, ...
    'SRmktAnnual',    SRmktAnnual, ...
    'SR',             SR, ...
    'SRAnnual',       SRAnnual, ...
    'appRatio',       appRatio, ...
    'appRatioAnnual', appRatioAnnual, ...
    'alphaAnnual',    alphaAnnual, ...
    'alphaAnnualVol', alphaAnnualVol, ...
    'alphaAnnualMom', alphaAnnualMom, ...
    'alphaAnnual3Fac', alphaAnnual3Fac, ...
    'alphaAnnualLPM',  alphaAnnualLPM);
end
