function n = useParfor(cfg)
% USEPARFOR  Return the parfor max-worker count, falling back to 0 (serial)
% when the Parallel Computing Toolbox is unavailable.
%
%   n = useParfor()        autodetect (recommended)
%   n = useParfor(cfg)     respect cfg.useParallel
%
% Use as the second argument to parfor:
%
%   parfor (k = 1:N, useParfor())
%       ...
%   end
%
% When n = 0 the parfor degrades to an ordinary for. This is the
% supported pattern for code that should run identically on machines
% without the Parallel Computing Toolbox.

if nargin >= 1 && isstruct(cfg) && isfield(cfg, 'useParallel') ...
        && ~isempty(cfg.useParallel)
    if cfg.useParallel
        n = Inf;
    else
        n = 0;
    end
    return
end

% Auto-detect
hasLicense = license('test', 'Distrib_Computing_Toolbox');
hasInstall = ~isempty(ver('parallel'));
if hasLicense && hasInstall
    n = Inf;
else
    n = 0;
end
end
