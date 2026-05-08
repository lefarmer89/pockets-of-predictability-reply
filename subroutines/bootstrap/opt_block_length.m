function Bstar = opt_block_length(data)
% function Bstar = opt_block_length(data);
%
% This is a function to select the optimal (in the sense
% of minimising the MSE of the estimator of the long-run
% variance) block length for the stationary bootstrap or circular bootstrap.
% Code follows Politis and White, 2001, "Automatic Block-Length
% Selection for the Dependent Bootstrap".
%
%  DECEMBER 2007: CORRECTED TO DEAL WITH ERROR IN LAHIRI'S PAPER, PUBLISHED
%  BY NORDMAN IN THE ANNALS OF STATISTICS
%
%  NOTE:    The optimal average block length for the stationary bootstrap, and it does not need to be an integer. 
%           The optimal block length for the circular bootstrap should be an integer. Politis and White suggest 
%               rounding the output UP to the nearest integer.
%
% INPUTS:	data, an nxk matrix
%
% OUTPUTS:	Bstar, a 2xk vector of optimal bootstrap block lengths, [BstarSB;BstarCB]
%
%  Andrew Patton
%
%  4 December, 2002
% 	Revised (to include CB): 13 January, 2003.
%
% Helpful suggestions for this code were received from Dimitris Politis and Kevin Sheppard.
%
%Modified 23.8.2003 by Kevin Sheppard for speed issues

[n,k] = size(data);

% these are optional in opt_block_length_full.m, but fixed at default values here
KN=max(5,sqrt(log10(n)));
%mmax = ceil(sqrt(n));
mmax = ceil(sqrt(n))+KN;           % adding KN extra lags to employ Politis' (2002) suggestion for finding largest signif m
warning_flags=0;
round=0;
%Bmax = sqrt(n);                  % maximum value of Bstar to consider.
Bmax = ceil(min(3*sqrt(n),n/3));  % dec07: new idea for rule-of-thumb to put upper bound on estimated optimal block length

c=2;
origdata=data;
Bstar_final=[];

for i=1:k
   data=origdata(:,i);
   
   % FIRST STEP: finding mhat-> the largest lag for which the autocorrelation is still significant.
   temp = mlag(data,mmax);
   temp = temp(mmax+1:end,:);	% dropping the first mmax rows, as they're filled with zeros
   temp = corrcoef([data(mmax+1:end),temp]);
   temp = temp(2:end,1);
   
   % We follow the empirical rule suggested in Politis, 2002, "Adaptive Bandwidth Choice".
   % as suggested in Remark 2.3, setting c=2, KN=5
   temp2 = [mlag(temp,KN)',temp(end-KN+1:end)];		% looking at vectors of autocorrels, from lag mhat to lag mhat+KN
   temp2 = temp2(:,KN+1:end);		% dropping the first KN-1, as the vectors have empty cells
   temp2 = (abs(temp2)<(c*sqrt(log10(n)/n)*ones(KN,mmax-KN+1)));	% checking which are less than the critical value
   temp2 = sum(temp2)';		% this counts the number of insignif autocorrels
   temp3 = [(1:1:length(temp2))',temp2];
   temp3 = temp3(find(temp2==KN),:);	% selecting all rows where ALL KN autocorrels are not signif
   if isempty(temp3)
      mhat = max(find(abs(temp) > (c*sqrt(log10(n)/n)) )); % this means that NO collection of KN autocorrels were all insignif, so pick largest significant lag
   else
      mhat = temp3(1,1);	% if more than one collection is possible, choose the smallest m
   end
   if 2*mhat>mmax;
      M = mmax;
      trunc1=1;
   else
      M = 2*mhat;
   end
   clear temp temp2 temp3;
   
   
   % SECOND STEP: computing the inputs to the function for Bstar
   kk = (-M:1:M)';
   
   if M>0;
      temp = mlag(data,M);
      temp = temp(M+1:end,:);	% dropping the first mmax rows, as they're filled with zeros
      temp = cov([data(M+1:end),temp]);
      acv = temp(:,1);			% autocovariances
      acv2 = [-(1:1:M)',acv(2:end)];
      if size(acv2,1)>1;
         acv2 = sortrows(acv2,1);
      end
      acv = [acv2(:,2);acv];			% autocovariances from -M to M
      clear acv2;
      Ghat = sum(lam(kk/M).*abs(kk).*acv);
      DCBhat = 4/3*sum(lam(kk/M).*acv)^2;
 
% OLD nov07      
%      DSBhat = 2/pi*quadl('opt_block_length_calc',-pi,pi,[],[],kk,acv,lam(kk/M));
%      DSBhat = DSBhat + 4*sum(lam(kk/M).*acv)^2;	% first part of DSBhat (note cos(0)=1)

% NEW dec07
      DSBhat = 2*(sum(lam(kk/M).*acv)^2);	% first part of DSBhat (note cos(0)=1)
      
      % FINAL STEP: constructing the optimal block length estimator
      Bstar = ((2*(Ghat^2)/DSBhat)^(1/3))*(n^(1/3));
      if Bstar>Bmax
         Bstar = Bmax;
      end
      BstarCB = ((2*(Ghat^2)/DCBhat)^(1/3))*(n^(1/3));
      
      if BstarCB>Bmax
         BstarCB = Bmax;
      end      
      Bstar = [Bstar;BstarCB];
   else
      Ghat = 0;
      % FINAL STEP: constructing the optimal block length estimator
      Bstar = [1;1];
   end
   Bstar_final=[Bstar_final Bstar];
end
Bstar=Bstar_final;

%%%%%%%%%%%%%%%%%%%%%%%%
function lam=lam(kk)
%Helper function, calculates the flattop kernel weights
lam = (abs(kk)>=0).*(abs(kk)<0.5)+2*(1-abs(kk)).*(abs(kk)>=0.5).*(abs(kk)<=1);

function xlag = mlag(x,n,init)
% PURPOSE: generates a matrix of n lags from a matrix (or vector)
%          containing a set of vectors (For use in var routines)
%---------------------------------------------------
% USAGE:     xlag = mlag(x,nlag)
%       or: xlag1 = mlag(x), which defaults to 1-lag
% where: x = a matrix (or vector), nobs x nvar
%     nlag = # of contiguous lags for each vector in x
%     init = (optional) scalar value to feed initial missing values
%            (default = 0)
%---------------------------------------------------
% RETURNS:
%        xlag = a matrix of lags (nobs x nvar*nlag)
%        x1(t-1), x1(t-2), ... x1(t-nlag), x2(t-1), ... x2(t-nlag) ...
% --------------------------------------------------
% SEE ALSO: lag() which works more conventionally
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

if nargin ==1 
n = 1; % default value
init = 0;
elseif nargin == 2
init = 0;
end;

if nargin > 3
error('mlag: Wrong # of input arguments');
end;

[nobs, nvar] = size(x);

xlag = ones(nobs,nvar*n)*init;
icnt = 0;
for i=1:nvar;
for j=1:n;
xlag(j+1:nobs,icnt+j) = x(1:nobs-j,i);
end;
icnt = icnt+n;
end;
