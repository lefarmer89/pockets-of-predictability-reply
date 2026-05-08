function [bsdata, indices] = stationary_bootstrap(data, B, w)
%STATIONARY_BOOTSTRAP  Politis & Romano (1994) stationary block bootstrap.
%
%   [BSDATA, INDICES] = stationary_bootstrap(DATA, B, W)
%
% Inputs:
%   DATA   T x 1 column vector to be bootstrapped.
%   B      Number of bootstrap samples (positive integer).
%   W      Average block length. Probability of starting a new block at
%          each step is P = 1/W.
%
% Outputs:
%   BSDATA   T x B matrix of bootstrapped values, BSDATA(:,j) = DATA(INDICES(:,j)).
%   INDICES  T x B index matrix into DATA (with wraparound). INDICES may
%            contain values up to 2*T-1; the caller's data array is
%            assumed to be doubled to accommodate this — the convention
%            from Kevin Sheppard's UCSD GARCH toolbox version, which the
%            existing dailyBootstrap caller already follows
%            (calib.residExtend = [resids; resids]).
%
% Replaces the per-iteration sequential loop with a fully vectorized
% cummax+broadcast computation. RNG sequence preserved (same calls in the
% same order: rand(1,B), rand(t,B), rand(1,nnz(select))) so the function
% is bit-identical to the original UCSD GARCH implementation for any
% fixed seed.
%
% Reference: Politis, D.N. and J.P. Romano (1994), "The stationary
% bootstrap." JASA 89(428), 1303-1313.

if nargin ~= 3
    error('stationary_bootstrap:TooFewInputs', '3 inputs required');
end

[t, k] = size(data);
if k > 1
    error('stationary_bootstrap:BadData', 'DATA must be a column vector');
end
if t < 2
    error('stationary_bootstrap:BadData', 'DATA must have at least 2 observations.');
end
if ~isscalar(w) || w < 1 || floor(w) ~= w
    error('stationary_bootstrap:BadW', 'W must be a positive scalar integer');
end
if ~isscalar(B) || B < 1 || floor(B) ~= B
    error('stationary_bootstrap:BadB', 'B must be a positive scalar integer');
end

p = 1 / w;

% Step 1: initial random positions for row 1.
indices = zeros(t, B);
indices(1, :) = ceil(t * rand(1, B));

% Step 2: per-cell block-restart selector (Bernoulli(p)).
select = rand(t, B) < p;

% Step 3: random restart values at every select=true position. Drawn in
% column-major order to match the original `indices(select) = ...` linear
% assignment.
indices(select) = ceil(rand(1, nnz(select)) * t);

% Step 4: vectorized fill of non-restart positions. Each column of
% `indices` is a piecewise sequence: at every block-start row, the value
% is a fresh random draw; between block-starts, it increments by 1.
%
% The original loop does this serially via
%   for i = 2:t
%       indices(i, ~select(i,:)) = indices(i-1, ~select(i,:)) + 1;
%   end
%
% Equivalent vectorized form: for each (i, j), find the row of the most
% recent block-start at or before i (call it `r`), then
%   indices(i, j) = indices(r, j) + (i - r).
% `r` is computed with cummax of (rowIdx .* sm), where sm marks block-start
% rows. Row 1 is treated as a block start (its random draw was set in
% step 1 / overwritten in step 3 if select(1,:) was true).
sm = select;
sm(1, :) = true;
rowIdx        = (1:t)';
blockStartRow = cummax(rowIdx .* sm, 1);
linIdx        = blockStartRow + (0:B-1) * t;
blockStartVal = indices(linIdx);
indices       = blockStartVal + (rowIdx - blockStartRow);

% Step 5: bsdata via mod-based wraparound (avoids materializing [data; data]
% which doubles memory for no gain). Bit-identical to the original
% `data = [data; data]; bsdata = data(indices)` because data is periodic
% modulo t.
if nargout >= 1
    bsdata = data(mod(indices - 1, t) + 1);
end
end
