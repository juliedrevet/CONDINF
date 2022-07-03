function [out] = fit_overall_stab_expe2(cfg)
%  FIT_STAB_EXPE2  Fit stabilization parameters fitting overall repetition
%  fraction from EXPE2 (=> 3 stimulus strength for simulation)
%
%  This function requires BADS v1.0 (https://github.com/lacerbi/bads/).

% check presence of required parameters
if ~all(isfield(cfg,{'t','s','r','w','pfal'}))
    error('Incomplete data information!');
end
% check presence of optional parameters
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
    fprintf('Assuming %d samples for the particle filter.\n',cfg.nsmp);
end
if ~isfield(cfg,'ftype')
    cfg.ftype = 'negllh';
    fprintf('Assuming %s as function value to be minimized.\n',cfg.ftype);
end
if ~isfield(cfg,'verbose')
    cfg.verbose = false;
end

% check function value type
if ~ismember(cfg.ftype,{'sse','negllh'})
    error('Undefined function value type: %s!',cfg.ftype);
end

% create data structure
dat   = [];
dat.t = cfg.t; % trial number in block
dat.s = cfg.s; % stimulus (+1/-1)
dat.r = cfg.r; % response (+1/-1)
dat.w = cfg.w; % stimulus weight
dat.r(isnan(dat.r)) = 0;

% check number of stimulus weights
if numel(unique(dat.w)) ~= 3
    error('Invalid number of stimulus weights!');
end

modtype = cfg.modtype;

ntrl    = numel(dat.t); % number of trials
nsmp    = cfg.nsmp; % number of samples
nval    = 100; % number of function evaluations at end point
ftype   = cfg.ftype; % function value type
verbose = cfg.verbose; % display level

% define useful function handles
logit = @(x)log(x./(1-x)); % logit function

% compute sensory noise
pfal = cfg.pfal; % false-alarm rate
if numel(pfal) ~= 3
    error('Invalid number of false-alarm rates!');
end
sigsen = 1/norminv(1-pfal(2)); % sensory noise
strsen = norminv(1-pfal)/norminv(1-pfal(2)); % sensory strengths
for iw = 1:3
    dat.s(dat.w == iw) = strsen(iw)*dat.s(dat.w == iw);
end

% compute participant repetition score
rep = [nan;dat.r(2:end) == dat.r(1:end-1)];
rep(dat.t == 1) = nan;
prep_sub = nanmean(rep);

% define model parameters
pnam = {}; % name
pmin = []; % minimum value
pmax = []; % maximum value
pfun = {}; % log-prior function
pini = []; % initial value
pplb = []; % plausible lower bound
ppub = []; % plausible upper bound
% define general parameters:
% 1/ perceived hazard rate
pnam{1,1} = 'h';
pmin(1,1) = 0;
pmax(1,1) = 0.5;
pfun{1,1} = @(x)betapdf(x,2,18);
pini(1,1) = betastat(2,18);
pplb(1,1) = betainv(0.15,2,18);
ppub(1,1) = betainv(0.85,2,18);
% 2/ inference noise
pnam{1,2} = 'siginf';
pmin(1,2) = 0.001;
pmax(1,2) = 5;
pfun{1,2} = @(x)gampdf(x,4,0.25);
pini(1,2) = gamstat(4,0.25);
pplb(1,2) = gaminv(0.15,4,0.25);
ppub(1,2) = gaminv(0.85,4,0.25);
% 3/ selection noise
pnam{1,3} = 'sigsel';
pmin(1,3) = 0.001;
pmax(1,3) = 5;
pfun{1,3} = @(x)gampdf(x,4,0.25);
pini(1,3) = gamstat(4,0.25);
pplb(1,3) = gaminv(0.15,4,0.25);
ppub(1,3) = gaminv(0.85,4,0.25);
% define model-specific parameters:
switch modtype
    case 'senbias' % sensory bias
        % 4/ sensory bias strength
        pnam{1,4} = 'lambda';
        pmin(1,4) = 0;
        pmax(1,4) = 5;
        pfun{1,4} = @(x)gampdf(x,1,0.75);
        pini(1,4) = gamstat(1,0.75);
        pplb(1,4) = gaminv(0.15,1,0.75);
        ppub(1,4) = gaminv(0.85,1,0.75);
    case 'inflaps' % inference lapses
        pnam{1,4} = 'plaps';
        pmin(1,4) = 0;
        pmax(1,4) = 0.999;
        pfun{1,4} = @(x)betapdf(x,1,9);
        pini(1,4) = betastat(1,9);
        pplb(1,4) = betainv(0.15,1,9);
        ppub(1,4) = betainv(0.85,1,9);
    case 'infdisc' % inference discard
        % 4/ discard criterion
        pnam{1,4} = 'delta';
        pmin(1,4) = 0;
        pmax(1,4) = 5;
        pfun{1,4} = @(x)gampdf(x,1,0.75);
        pini(1,4) = gamstat(1,0.75);
        pplb(1,4) = gaminv(0.15,1,0.75);
        ppub(1,4) = gaminv(0.85,1,0.75);
    case 'selbias' % selection bias
        % 4/ selection bias strength
        pnam{1,4} = 'beta';
        pmin(1,4) = 0;
        pmax(1,4) = 2;
        pfun{1,4} = @(x)gampdf(x,1,0.25);
        pini(1,4) = gamstat(1,0.25);
        pplb(1,4) = gaminv(0.15,1,0.25);
        ppub(1,4) = gaminv(0.85,1,0.25);
    case 'selepsi' % epsilon-greedy selection
        pnam{1,4} = 'epsi';
        pmin(1,4) = 0;
        pmax(1,4) = 0.999;
        pfun{1,4} = @(x)betapdf(x,1,9);
        pini(1,4) = betastat(1,9);
        pplb(1,4) = betainv(0.15,1,9);
        ppub(1,4) = betainv(0.85,1,9);
    case 'idindep' % independant inference discard
        % 4/ discard criterion 1
        pnam{1,4} = 'delta1';
        pmin(1,4) = 0;
        pmax(1,4) = 5;
        pfun{1,4} = @(x)gampdf(x,1,0.75);
        pini(1,4) = gamstat(1,0.75);
        pplb(1,4) = gaminv(0.15,1,0.75);
        ppub(1,4) = gaminv(0.85,1,0.75);
        % 5/ discard criterion 2
        pnam{1,5} = 'delta2';
        pmin(1,5) = 0;
        pmax(1,5) = 5;
        pfun{1,5} = @(x)gampdf(x,1,0.75);
        pini(1,5) = gamstat(1,0.75);
        pplb(1,5) = gaminv(0.15,1,0.75);
        ppub(1,5) = gaminv(0.85,1,0.75);
        % 6/ discard criterion 3
        pnam{1,6} = 'delta3';
        pmin(1,6) = 0;
        pmax(1,6) = 5;
        pfun{1,6} = @(x)gampdf(x,1,0.75);
        pini(1,6) = gamstat(1,0.75);
        pplb(1,6) = gaminv(0.15,1,0.75);
        ppub(1,6) = gaminv(0.85,1,0.75);
end

% check fixed parameters
% check general parameters:
if isfield(cfg,'h')
    cfg.h = min(max(cfg.h,0),0.5);
end
if isfield(cfg,'siginf')
    cfg.siginf = max(cfg.siginf,1e-6);
end
if isfield(cfg,'sigsel')
    cfg.sigsel = max(cfg.sigsel,1e-6);
end
% check model-specific parameters:
switch modtype
    case 'senbias'
        if isfield(cfg,'lambda')
            cfg.lambda = max(cfg.lambda,0);
        end
    case 'inflaps'
        if isfield(cfg,'plaps')
            cfg.plaps = min(max(cfg.plaps,0),1);
        end
    case 'infdisc'
        if isfield(cfg,'delta')
            cfg.delta = max(cfg.delta,0);
        end
    case 'selbias'
        if isfield(cfg,'beta')
            cfg.beta = max(cfg.beta,0);
        end
    case 'selepsi'
        if isfield(cfg,'epsi')
            cfg.epsi = min(max(cfg.epsi,0),1);
        end
    case 'idindep'
        if isfield(cfg,'delta1')
            cfg.delta1 = max(cfg.delta1,0);
        end
        if isfield(cfg,'delta2')
            cfg.delta2 = max(cfg.delta2,0);
        end
        if isfield(cfg,'delta3')
            cfg.delta3 = max(cfg.delta3,0);
        end
end

% define fixed parameters
npar = numel(pnam);
pfix = cell(1,npar);
for i = 1:npar
    if isfield(cfg,pnam{i})
        pfix{i} = cfg.(pnam{i});
    end
end

% define free parameters
ifit = cell(1,npar);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
pfit_plb = [];
pfit_pub = [];
n = 1;
for ipar = 1:npar
    if isempty(pfix{ipar}) % free parameter
        ifit{ipar} = n;
        pfit_ini = cat(2,pfit_ini,pini(ipar));
        pfit_min = cat(2,pfit_min,pmin(ipar));
        pfit_max = cat(2,pfit_max,pmax(ipar));
        pfit_plb = cat(2,pfit_plb,pplb(ipar));
        pfit_pub = cat(2,pfit_pub,ppub(ipar));
        n = n+1;
    end
end
nfit = length(pfit_ini);

if nfit > 0
    
    % set BADS options
    options = bads('defaults');
    options.UncertaintyHandling = 1; % noisy objective function
    options.NoiseFinalSamples = nval; % number of samples
    % set display level
    if verbose
        options.Display = 'iter';
    else
        options.Display = 'none';
    end
    
    % fit model using BADS
    [phat,fval,exitflag,output] = bads(@(x)fun(x), ...
        pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,[],options);
    
    % create output structure with parameter values
    [~,phat] = fun(phat);
    out = cell2struct(phat(:),pnam(:));
    
    out.modtype = modtype;
    
    % store function value information
    out.ftype = ftype; % function value type
    out.fval  = output.fval; % estimated function value
    out.fstd  = output.fsd; % estimated s.d. of function value
    
    % store additional output from BADS
    out.exitflag = exitflag;
    out.output = output;
    
else
    
    % create output structure with parameter values
    [~,phat] = fun([]);
    out = cell2struct(phat(:),pnam(:));
    
    % estimate function value
    fval = nan(nval,1);
    for ival = 1:nval
        fval(ival) = fun([]);
    end
    
    % store function value information
    out.ftype = ftype; % function value type
    out.fval = mean(fval); % estimated function value
    out.fstd = std(fval); % estimated s.d. of function value
    
end

    function [f,pval] = fun(p)
        % retrieve parameter values
        pval = cell(1,npar);
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                pval{k} = p(ifit{k});
            else % fixed parameter
                pval{k} = pfix{k};
            end
        end
        % simulate reversal metrics
        resp = simr((pval{:}));
        
        rep_sim = [nan(1,nsmp);resp(2:end,:) == resp(1:end-1,:)];
        rep_sim(dat.t == 1,:) = nan;
        
        prep_avg = nanmean(rep_sim);
        prep_std = nanstd(prep_avg);
        % compute function value to be minimized
        switch ftype
            case 'sse'
                % compute sum of squared errors
                f = (rep_avg-prep_sub).^2;
            case 'negllh'
                % compute negative log-likelihood
                seps = 0.001; % minimum standard deviation
                f = -normllh(prep_sub,mean(prep_avg),max(prep_std,seps));
        end
    end

   function [r] = simr(varargin)
        switch modtype
            case 'flatinf' % flat inference
                r = simr_flatinf(varargin{:});
            case 'senbias' % sensory bias
                r = simr_senbias(varargin{:});
            case 'inflaps' % inference lapses
                r = simr_inflaps(varargin{:});
            case 'infdisc' % inference discard
                r = simr_infdisc(varargin{:});
            case 'selbias' % selection bias
                r = simr_selbias(varargin{:});
            case 'selepsi' % epsilon-greedy selection
                r = simr_selepsi(varargin{:});
            case 'idindep' % inndependant discard criteria
                r = simr_idindep(varargin{:});
        end
    end

    function [r] = simr_flatinf(h,siginf,sigsel)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % add sensory stimulus
            s = dat.s(itrl)+sigsen*randn(1,nsmp);
            s = s*2/sigsen^2;
            x(itrl,:) = x(itrl,:)+s;
            % add inference noise
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
    end

    function [r] = simr_senbias(h,siginf,sigsel,lambda)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = dat.s(itrl)+sigsen*randn(1,nsmp);
            s = s*2/sigsen^2;
            ip = s+x(itrl,:)*lambda > 0;
            in = ~ip;
            % add log-evidence
            x(itrl,ip) = x(itrl,ip)+ ...
                log(1-normcdf(x(itrl,ip)*lambda,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(x(itrl,ip)*lambda,-2/sigsen^2,2/sigsen));
            x(itrl,in) = x(itrl,in)+ ...
                log(normcdf(x(itrl,in)*lambda,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(x(itrl,in)*lambda,-2/sigsen^2,2/sigsen));
            % add inference noise
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
    end

    function [r] = simr_inflaps(h,siginf,sigsel,plaps)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = sign(dat.s(itrl)+sigsen*randn(1,nsmp));
            % sample inference lapses
            i0 = rand(1,nsmp) < plaps;
            % add log-evidence for non-lapsing trials
            x(itrl,~i0) = x(itrl,~i0)+s(~i0)*logit(1-pfal(2));
            % add inference noise for non-lapsing trials
            x(itrl,~i0) = x(itrl,~i0)+siginf*randn(1,nnz(~i0));
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
    end

    function [r] = simr_infdisc(h,siginf,sigsel,delta)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = dat.s(itrl)+sigsen*randn(1,nsmp);
            s = s*2/sigsen^2;
            ip = s > +delta;
            in = s < -delta;
            i0 = ~ip & ~in;
            % add log-evidence
            x(itrl,ip) = x(itrl,ip)+ ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen));
            x(itrl,in) = x(itrl,in)+ ...
                log(normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add inference noise
            x(itrl,~i0) = x(itrl,~i0)+siginf*randn(1,nnz(~i0));
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
    end

    function [r] = simr_selbias(h,siginf,sigsel,beta)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = sign(dat.s(itrl)+sigsen*randn(1,nsmp));
            % add log-evidence
            x(itrl,:) = x(itrl,:)+s*logit(1-pfal(2));
            % add inference noise
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response with selection noise
            if dat.t(itrl) == 1
                % w/o selection bias
                r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
            else
                % with selection bias
                r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp)+beta*r(itrl-1,:));
            end
        end
    end

    function [r] = simr_selepsi(h,siginf,sigsel,epsi)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        for itrl = 1:ntrl
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = sign(dat.s(itrl)+sigsen*randn(1,nsmp));
            % add log-evidence
            x(itrl,:) = x(itrl,:)+s*logit(1-pfal(2));
            % add inference noise
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response
            if dat.t(itrl) == 1
                % sample response with selection noise
                r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
            else
                i0 = rand(1,nsmp) < epsi;
                % blind response repetition
                r(itrl,i0) = r(itrl-1,i0);
                % sample response with selection noise
                r(itrl,~i0) = sign(x(itrl,~i0)+sigsel*randn(1,nnz(~i0)));
            end
        end
    end

    function [r] = simr_idindep(h,siginf,sigsel,delta1,delta2,delta3)
        x = zeros(ntrl,nsmp); % simulated log-belief
        r = zeros(ntrl,nsmp); % simulated response
        delta_vec = [delta1 delta2 delta3];
        for itrl = 1:ntrl
            delta = delta_vec(dat.w(itrl)); % trial specific criterion
            % add log-prior
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
            end
            % categorize sensory stimulus
            s = dat.s(itrl)+sigsen*randn(1,nsmp);
            s = s*2/sigsen^2;
            ip = s > +delta;
            in = s < -delta;
            i0 = ~ip & ~in;
            % add log-evidence
            x(itrl,ip) = x(itrl,ip)+ ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen));
            x(itrl,in) = x(itrl,in)+ ...
                log(normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add inference noise
            x(itrl,~i0) = x(itrl,~i0)+siginf*randn(1,nnz(~i0));
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
    end

end

function [y] = update_prior(x,h)
% update prior belief wrt hazard rate
y = x+log((1-h)./h+exp(-x))-log((1-h)./h+exp(+x));
end
