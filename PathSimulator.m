function out = PathSimulator(theta, TSample, random_seed, lambda, outMatfile)

% load calibration
thetaBest = theta;

%% SET UP PARAMETERS
N = 7; % number of state variables
K = 1;  % number of states
e = eye(N); % basis vectors
mu_rf = .015 / 252; % mean risk-free rate
rhoCF = thetaBest(9); %0.469^(1/252);  % cash-flow state persistence
mu_r  = .085 / 252; % mean returns
mu_rp = mu_r - mu_rf; % mean risk-premium
mean_dp = -3.43 - log(252); % mean log dividend-price ratio

% linearization constants
rho = 1/(1+exp(mean_dp));
k = -log(rho) - (1-rho)*log(1/rho-1);

% mean of dividend growth
mu_d = mu_r - (k - rho * mean_dp + mean_dp);

%% CALIBRATED PARAMETERS

rhoDR      = thetaBest(1); % 1. persistence of discount rate state
rhoTP      = thetaBest(2); % 2. persistence of time pref state
beta_rf_cf = thetaBest(3); % 3. loading of risk-free on cf state
beta_rf_dr = thetaBest(4); % 4. loading of risk-free on dr state
sigma_cf   = thetaBest(5); % 5. volatility of shock to cf
sigma_dr   = thetaBest(6); % 6. volatility of shock to dr
sigma_tp   = thetaBest(7); % 7. volatility of shock to tp
sigma_dg   = thetaBest(8); % 8. volatility of shock to dividend growth

%% MEAN VECTOR

% preallocate
muYMat = zeros(N,K);

% mean of expected dividend growth
muYMat(4,:) = (1-lambda) * mu_d;

% mean of dividend growth
muYMat(5,:) = mu_d;

% mean of risk premium
muYMat(6,:) = mu_rp;

% mean of risk-free rate
muYMat(7,:) = mu_rf - lambda * beta_rf_cf * mu_d;

%% AR COEFFICIENT MATRIX

% preallocate
AMat = zeros(N,N,K);

% persistence of state variables
AMat(1:3,1:3,:) = repmat(diag([rhoCF, rhoDR, rhoTP]),[1 1 K]);

% expected dividend growth
AMat(4,1,:) = (1-lambda) * rhoCF;
AMat(4,4,:) = lambda * rhoCF;

% dividend growth
AMat(5,1,:) = 1;

% expected risk premium
AMat(6,2,:) = rhoDR;

% risk-free rate
AMat(7,1,:) = (1 - lambda) * rhoCF * beta_rf_cf;
AMat(7,2,:) = rhoDR * beta_rf_dr;
AMat(7,3,:) = rhoTP;
AMat(7,4,:) = lambda * rhoCF * beta_rf_cf;

%% VARIANCE-COVARIANCE MATRIX OF SHOCKS

% preallocate
SigMat = zeros(N,N,K);

% shocks to state variables
SigMat(1,1,:) = sigma_cf;
SigMat(2,2,:) = sigma_dr;
SigMat(3,3,:) = sigma_tp;

% shock to expected dividend growth
SigMat(4,:,:) = (1-lambda) * SigMat(1,:,:);

% shock to realized dividend growth
SigMat(5,5,:) = sigma_dg;

% shock to expected risk premium
SigMat(6,:,:) = SigMat(2,:,:);

% shock to risk-free rate
for k = 1:K
    SigMat(7,:,k) = ((1-lambda) * beta_rf_cf * e(1,:) + ...
        beta_rf_dr * e(2,:) + e(3,:)) * SigMat(:,:,k);
end

%% Analytic variances

% useful selector matrices
Spd = [0; -(1 + beta_rf_dr) / (1 - rho * rhoDR); ...
    -1 / (1 - rho * rhoTP); (1 - beta_rf_cf) / (1-rho*rhoCF); ...
    zeros(3,1)];

%% Main loop

% TSample = 23786;
T = TSample; %5000000; % number of daily observations to simulate

% random seed
rng(random_seed)

% states
YSim = [muYMat(:,1),NaN(N,T)];

% linearization constants
rho = 1/(1+exp(mean_dp));
k = -log(rho) - (1-rho)*log(1/rho-1);

dpSim = NaN(T+1,1);  % log dividend-price ratio
rSim = [mu_r; NaN(T,1)];

eSim = (SigMat * mvnrnd(zeros(1,N), diag(ones(1,N)), T+1)');

for t = 1:(T+1)

    if ~mod(t,TSample*5)
        round(t/TSample);
    end

    % Bianchi expectations

    dpSim(t) = -k/(1-rho) - (mu_d - mu_rf - mu_rp)/(1-rho) + ...
        mu_d * (1-beta_rf_cf)/(1-rho*rhoCF) - Spd'*YSim(:,t);

    if t > 1 && t < T+1
        rSim(t) = k - rho*dpSim(t) + YSim(N-2,t) + dpSim(t-1);
    end

    % Simulate next period

    if t < T+1

        YSim(:,t+1) = muYMat + AMat*YSim(:,t) + eSim(:,t+1);

    end

end

% Stagger variables appropriately to correct for t+2 timing of rf
dpSim = dpSim(2:end);
dSim = YSim(N-2,2:end)';
rfSim = YSim(N,1:end-1)';
rpSim = YSim(N-1,2:end)';
rSim = rSim(2:end);

% outcome variable (excess returns)
% if minus_ZDR == 1
%     erSim = rSim-rfSim - YSim(2,2:end)';
% else
erSim = rSim-rfSim;
% end

zSim = YSim(1:(N-3),2:end)';
sSim = ones(size(dpSim));

% construct rvar predictor
rvarSim = rSim.^2;
squaredret = rSim.^2;
for tt = 61:length(rvarSim)
    rvarSim(tt) =sum(squaredret((tt-60):tt), 'omitnan');
end

save(outMatfile, 'YSim' ,'dpSim' , 'dSim', 'rfSim', 'rpSim', 'rSim',...
        'erSim', 'zSim', 'sSim', 'rvarSim');

out = YSim;

end