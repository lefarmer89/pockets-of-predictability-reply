function [fCW, fDM] = forecastDifferentials(yActual, yF1, yF2)
% FORECASTDIFFERENTIALS  Per-period Clark-West and Diebold-Mariano
% forecast differentials.
%
%   [fCW, fDM] = forecastDifferentials(yActual, yF1, yF2)
%
% Inputs:
%   yActual  T x 1 realized values
%   yF1      T x 1 benchmark forecast
%   yF2      T x 1 challenger forecast
%
% Outputs:
%   fCW      Clark-West differential
%   fDM      Diebold-Mariano differential
%
% Construction matches the convention used throughout the empirics:
%   fDM = (y - yF1)^2 - (y - yF2)^2
%   fCW = (y - yF1)^2 - ((y - yF2)^2 - (yF1 - yF2)^2)

fDM = (yActual - yF1).^2 - (yActual - yF2).^2;
fCW = (yActual - yF1).^2 - ((yActual - yF2).^2 - (yF1 - yF2).^2);
end
