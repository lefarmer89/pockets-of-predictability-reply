function [K, k2] = oneSidedKernel()
% ONESIDEDKERNEL  The one-sided Epanechnikov kernel used throughout the
% pipeline, plus its second-moment constant k2.
%
%   [K, k2] = oneSidedKernel()
%
%   K(z)  = 1.5 * (1 - z^2) for -1 <= z <= 0, else 0
%   k2    second-moment constant (1.2 for the one-sided Epanechnikov)
%
% Centralizing this lets every entry-point use the same kernel without
% re-defining it inline.

K = @(z) 1.5 .* (1 - z.^2) .* (abs(z) <= 1) .* (z <= 0);
k2 = 1.2;
end
