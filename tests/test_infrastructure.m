function r = test_infrastructure()
% TEST_INFRASTRUCTURE  Verify setup_paths, default_config, and useParfor
% behave correctly. These three helpers are the foundation for cfg-driven
% entry-point functions in Phase 6+.

r.name = 'test_infrastructure';
r.pass = false;
r.message = '';

checks = {};

% ---- setup_paths ----
try
    paths = setup_paths();
    requiredFields = {'root','subroutines','tables','figures','data', ...
        'csvSims','results','output','outputTables','outputFigures', ...
        'outputLogs'};
    for k = 1:numel(requiredFields)
        if ~isfield(paths, requiredFields{k})
            r.message = sprintf('paths.%s missing', requiredFields{k});
            return
        end
    end
    if ~isfolder(paths.root)
        r.message = sprintf('paths.root does not exist: %s', paths.root);
        return
    end
    checks{end+1} = sprintf('setup_paths OK (%d fields)', numel(requiredFields));
catch ME
    r.message = sprintf('setup_paths failed: %s', ME.message);
    return
end

% ---- default_config ----
try
    cfg = default_config();
    requiredCfg = {'paths','coefWindowYears','sedWindowYears', ...
        'weightWindowYears','minPocketDays','gPrior','shrinkageTarget', ...
        'winsorPct','weightLB','weightUB','adjCostBps','sampleSplit', ...
        'selectionCriterion','useParallel','rngSeed'};
    for k = 1:numel(requiredCfg)
        if ~isfield(cfg, requiredCfg{k})
            r.message = sprintf('cfg.%s missing', requiredCfg{k});
            return
        end
    end
    % Verify defaults match the published reply
    if cfg.gPrior ~= 2
        r.message = sprintf('cfg.gPrior default should be 2, got %g', cfg.gPrior);
        return
    end
    if ~strcmp(cfg.shrinkageTarget, 'zero')
        r.message = sprintf('cfg.shrinkageTarget default should be ''zero''');
        return
    end
    if cfg.minPocketDays ~= 21
        r.message = sprintf('cfg.minPocketDays default should be 21');
        return
    end
    checks{end+1} = sprintf('default_config OK (%d fields)', numel(requiredCfg));
catch ME
    r.message = sprintf('default_config failed: %s', ME.message);
    return
end

% ---- useParfor ----
try
    n = useParfor();
    if ~(isnumeric(n) && isscalar(n) && (n == 0 || n == Inf))
        r.message = 'useParfor() should return 0 or Inf';
        return
    end
    n0 = useParfor(struct('useParallel', false));
    if n0 ~= 0
        r.message = 'useParfor(cfg.useParallel=false) should return 0';
        return
    end
    nInf = useParfor(struct('useParallel', true));
    if nInf ~= Inf
        r.message = 'useParfor(cfg.useParallel=true) should return Inf';
        return
    end
    checks{end+1} = sprintf('useParfor OK (auto=%g)', n);
catch ME
    r.message = sprintf('useParfor failed: %s', ME.message);
    return
end

% ---- entry-point function parse check ----
% Verify each converted entry point loads as a function (i.e. parses
% without syntax errors and reports its arity).
try
    pkgRoot = fileparts(fileparts(mfilename('fullpath')));
    entryPoints = {
        'dailyEmpirics', ...
        'dailyEmpiricsHyperparameters1', ...
        'dailyEmpiricsHyperparameters2', ...
        'dailyEmpiricsHyperparametersMarginals', ...
        'monthlyEmpirics', ...
        'famaFrenchEmpirics', ...
        'dailyBootstrap', ...
        'dailyAssetPricingBootstrap', ...
        'generateStickyExpectationsData', ...
        'stickyExpectationsSim'};
    nConverted = 0;
    nScript = 0;
    issues = {};
    for k = 1:numel(entryPoints)
        name = entryPoints{k};
        if exist(name, 'file') ~= 2
            issues{end+1} = [name ': not on path']; %#ok<AGROW>
            continue
        end
        try
            n = nargin(name);
            % nargin returns -1 for varargin; >= 0 means a real function
            if n >= 0
                nConverted = nConverted + 1;
            else
                nScript = nScript + 1;
            end
        catch parseErr
            % nargin throws if called on a script
            if contains(parseErr.message, 'is a script')
                nScript = nScript + 1;
            else
                issues{end+1} = sprintf('%s: %s', name, parseErr.message); %#ok<AGROW>
            end
        end
    end
    if isempty(issues)
        checks{end+1} = sprintf('entry-points: %d functions, %d scripts', ...
            nConverted, nScript);
    else
        r.message = strjoin(issues, '; ');
        return
    end
catch ME
    r.message = sprintf('entry-point parse check failed: %s', ME.message);
    return
end

r.pass = true;
r.message = strjoin(checks, '; ');
end
