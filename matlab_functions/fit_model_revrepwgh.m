function [out] = fit_model_revrepwgh(cfg)
%  FIT_MODEL_REVREPWGH  Fit model to reversal and repetition curves with
%  respect to stimulus weight
%
%  This function fits a desired model to reversal and repetition curves:
%    * modtype = 'flatinf' => flat inference
%    * modtype = 'nonstab' => hierarchical inference w/o stabilization
%    * modtype = 'senbias' => perceptual bias at sensory stage
%    * modtype = 'inflaps' => inference lapses
%    * modtype = 'infdisc' => conditional inference
%    * modtype = 'selbias' => repetition bias at selection stage
%    * modtype = 'selepsi' => response lapses : epsilon-greedy selection
%    * modtype = 'infwght' => weighted inference (added after review)
%
%  This function requires VBMC v1.0 (https://github.com/lacerbi/vbmc/).
%
%  Valentin Wyart <valentin.wyart@ens.fr> - June 2020 / Julie Drevet 2021

% check presence of required parameters
if ~all(isfield(cfg,{'t','s','w'}))
    error('Incomplete experiment information!');
end
if ~isfield(cfg,'modtype')
    error('Missing model type!');
end
if ~isfield(cfg,'pfal')
    error('Missing false-alarm rate!');
end
% check presence of optional parameters
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
end
if ~isfield(cfg,'std_sim')
    cfg.std_sim = true;
end
if ~isfield(cfg,'nres')
    cfg.nres = 1e3;
end
if ~isfield(cfg,'nval')
    cfg.nval = 0;
end
if ~isfield(cfg,'verbose')
    cfg.verbose = false;
end

% define useful function handles
logit = @(x)log(x./(1-x)); % logit function

% create data structure
dat = [];
dat.t = cfg.t; % trial number in block
dat.s = cfg.s; % stimulus: will be normalized according to strength
dat.c = cfg.s; % stimulus: as underlying category (+1/-1)
dat.w = cfg.w; % stimulus weight or strength (1:low|2:middle|3:high)

% check number of stimulus weights
if numel(unique(dat.w)) ~= 3
    error('Invalid number of stimulus weights!');
end

if all(isfield(cfg,{'crev_w','crep_w'}))
    if ...
            ndims(cfg.crev_w) ~= 2 || any(size(cfg.crev_w) ~= [8,3]) || ...
            ndims(cfg.crep_w) ~= 2 || any(size(cfg.crep_w) ~= [8,3])
        error('Invalid precomputed reversal/repetition curves!');
    end
    % use precomputed reversal/repetition curves
    use_cdat = true;
    dat.crev_w = cfg.crev_w; % reversal curve
    dat.crep_w = cfg.crep_w; % repetition curve
else
    % compute reversal/repetition curves from responses
    use_cdat = false;
    if ~isfield(cfg,'r')
        error('Missing response responses!');
    end
    dat.r = cfg.r; % response (+1/-1)
    dat.r(isnan(dat.r)) = 0; % time-out responses
end

% get number of reversals
nrev = nnz(dat.c(2:end) ~= dat.c(1:end-1) & dat.t(2:end) > 1);

% get model type
modtype = lower(cfg.modtype);
if ~ismember(modtype,{'flatinf','nonstab','senbias','inflaps','infdisc','selbias','selepsi','idindep','infwght'})
    error('Undefined model type!');
end
% fit non-stabilized inference as inference discarding w/o discarding
if strcmp(modtype,'nonstab')
    modtype = 'infdisc';
    cfg.delta = 0;
end

ntrl = numel(dat.t); % number of trials
std_sim = cfg.std_sim; % use simulations uncertainty?
nsmp = cfg.nsmp; % number of samples
nres = cfg.nres; % number of re-samples for bootstrapping
nval = cfg.nval; % number of function evaluations for validation
verbose = cfg.verbose; % display level

% compute sensory noise
pfal = cfg.pfal; % false-alarm rate
if numel(pfal) ~= 3
    error('Invalid number of false-alarm rates!');
end
sigsen = 1/norminv(1-pfal(2)); % sensory noise
strsen = norminv(1-pfal)/norminv(1-pfal(2)); % sensory strengths
% normalize stimulus
for iw = 1:3
    dat.s(dat.w == iw) = strsen(iw)*dat.s(dat.w == iw);
end

if use_cdat
    % use precomputed reversal/repetition curves
    crev_sub = dat.crev_w(:);
    crep_sub = dat.crep_w(:);
else
    % compute reversal/repetition curves
    csub = getc(dat.t,dat.c,dat.r,dat.w);
    crev_sub = csub.rev_w_avg(:); % reversal curve
    crep_sub = csub.rep_w_avg(:); % repetition curve
end

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
    case 'infwght' % weighted inference
        % 4/ weight criterion
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
    case 'infwght'
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
    
    % configure VBMC
    options = vbmc('defaults');
    options.MaxIter = 300; % maximum number of iterations
    options.MaxFunEvals = 500; % maximum number of function evaluations
    options.SpecifyTargetNoise = true; % noisy log-posterior function
    % set display level
    if verbose
        options.Display = 'iter';
    else
        options.Display = 'final';
    end
    
    % fit model using VBMC
    [vp,elbo,~,exitflag,output] = vbmc(@(x)fun(x), ...
        pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,options);
    
    % generate 10^6 samples from the variational posterior
    xsmp = vbmc_rnd(vp,1e6);
    
    xmap = vbmc_mode(vp); % posterior mode
    xavg = mean(xsmp,1); % posterior mean
    xstd = std(xsmp,[],1); % posterior s.d.
    xcov = cov(xsmp); % posterior covariance matrix
    xmed = median(xsmp,1); % posterior medians
    xiqr = quantile(xsmp,[0.25,0.75],1); % posterior interquartile ranges
    
    % create output structure with parameter values
    xhat = xavg; % use posterior average
    phat = getpval(xhat);
    out = cell2struct(phat(:),pnam(:));
    
    % create substructures with parameter values
    phat_map = getpval(xmap); % posterior mode
    out.pmap = cell2struct(phat_map(:),pnam(:));
    phat_avg = getpval(xavg); % posterior mean
    out.pavg = cell2struct(phat_avg(:),pnam(:));
    
    % store fitting information
    out.nsmp = nsmp; % number of samples used by particle filter
    out.nres = nres; % number of re-samples for bootstrapping
    
    % store VBMC output
    out.vp   = vp;
    out.elbo = elbo; % ELBO (expected lower bound on log-marginal likelihood)
    out.xnam = pnam(cellfun(@isempty,pfix)); % fitted parameters
    out.xmap = xmap; % posterior mode
    out.xavg = xavg; % posterior means
    out.xstd = xstd; % posterior s.d.
    out.xcov = xcov; % posterior covariance matrix
    out.xmed = xmed; % posterior medians
    out.xiqr = xiqr; % posterior interquartile ranges
    
    % store extra VBMC output
    out.exitflag = exitflag;
    out.output   = output;
    
else
    
    % create output structure with fixed parameter values
    phat = getpval([]);
    out = cell2struct(phat(:),pnam(:));
    
end

% remove useless parameter
if strcmpi(cfg.modtype,'nonstab')
    cfg = rmfield(cfg,'delta');
    out = rmfield(out,'delta');
end

% store subject reversal metrics
out.crev_w_sub = reshape(crev_sub,[8,3]); % reversal curve
out.crep_w_sub = reshape(crep_sub,[8,3]); % repetition curve

% store best-fitting reversal metrics
chat = simc(phat{:});
% reversal curve
out.crev_w_avg = chat.rev_w_avg;
out.crev_w_std = chat.rev_w_std;
% repetition curve
out.crep_w_avg = chat.rep_w_avg;
out.crep_w_std = chat.rep_w_std;

if nfit > 0 && nval > 0
    % perform function evaluations for validation
    l_val = nan(nval,1); % estimates for each validation sample
    l_sd_val = nan(nval,1); % bootstrap s.d. estimates for each validation sample
    for ival = 1:nval
        [l_val(ival),l_sd_val(ival)] = fun(xhat); % use posterior average
    end
    out.l_val = l_val; % estimates
    out.l_sd_val = l_sd_val; % bootstrap s.d. estimates
end

% store configuration structure
out.cfg = cfg;

    function [l,l_sd] = fun(p)
        % get parameter values
        pval = getpval(p);
        % get log-prior
        lp = 0;
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                lp = lp+log(pfun{k}(pval{k}));
            end
        end
        % get log-likelihood
        [ll,ll_sd] = getll(pval{:});
        % get log-posterior statistics
        l = sum(ll)+lp; % estimate
        l_sd = ll_sd; % bootstrap s.d. estimate
    end

    function [pval] = getpval(p)
        pval = cell(1,npar); % parameter values
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                pval{k} = p(ifit{k});
            else % fixed parameter
                pval{k} = pfix{k};
            end
        end
    end

    function [ll,ll_sd] = getll(varargin)
        % simulate reversal metrics
        c = simc(varargin{:});
        c.rev_w_all = reshape(c.rev_w_all,[],nsmp);
        c.rev_w_avg = c.rev_w_avg(:);
        c.rev_w_std = c.rev_w_std(:);
        c.rep_w_all = reshape(c.rep_w_all,[],nsmp);
        c.rep_w_avg = c.rep_w_avg(:);
        c.rep_w_std = c.rep_w_std(:);
        % compute bootstrap estimate of log-likelihood s.d.
        lres = nan(nres,1);
        for ires = 1:nres
            jres = randsample(nsmp,nsmp,true);
            crev_avg = mean(c.rev_w_all(:,jres),2);
            crep_avg = mean(c.rep_w_all(:,jres),2);
            if std_sim
                % use simulations uncertainty
                crev_std = std(c.rev_w_all(:,jres),[],2);
                crep_std = std(c.rep_w_all(:,jres),[],2);
            else
                % use data uncertainty
                crev_std = sqrt(crev_sub.*(1-crev_sub).*nrev);
                crep_std = sqrt(crep_sub.*(1-crep_sub).*nrev);
            end
            lres(ires) = ...
                sum(normllh(crev_sub,crev_avg,max(crev_std,1e-3)))+ ...
                sum(normllh(crep_sub,crep_avg,max(crep_std,1e-3)));
        end
        ll_sd = max(std(lres),1e-6);
        % compute log-likelihood estimate
        ll = ...
            sum(normllh(crev_sub,c.rev_w_avg,max(c.rev_w_std,1e-3)))+ ...
            sum(normllh(crep_sub,c.rep_w_avg,max(c.rep_w_std,1e-3)));
    end

    function [c] = simc(varargin)
        switch modtype
            case 'flatinf' % flat inference
                c = simc_flatinf(varargin{:});
            case 'senbias' % sensory bias
                c = simc_senbias(varargin{:});
            case 'inflaps' % inference lapses
                c = simc_inflaps(varargin{:});
            case 'infdisc' % inference discard
                c = simc_infdisc(varargin{:});
            case 'infwght' % weighted inference
                c = simc_infwght(varargin{:});
            case 'selbias' % selection bias
                c = simc_selbias(varargin{:});
            case 'selepsi' % epsilon-greedy selection
                c = simc_selepsi(varargin{:});
            case 'idindep' % inference discard
                c = simc_idindep(varargin{:});
        end
    end

    function [c] = simc_flatinf(h,siginf,sigsel)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_senbias(h,siginf,sigsel,lambda)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_inflaps(h,siginf,sigsel,plaps)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_infdisc(h,siginf,sigsel,delta)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_infwght(h,siginf,sigsel,delta)
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
            iphi = s > +delta; % above threshold
            inhi = s < -delta;
            iplo = s < +delta & s > 0; % below threshold
            inlo = s > -delta & s < 0;
            % add log-evidence => ABOVE THRESHOLD (hi)
            x(itrl,iphi) = x(itrl,iphi)+ ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen));
            x(itrl,inhi) = x(itrl,inhi)+ ...
                log(normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add log-evidence => BELOW THRESHOLD (lo)
            x(itrl,iplo) = x(itrl,iplo)+ ...
                log(normcdf(+delta,+2/sigsen^2,2/sigsen)-normcdf(0,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(+delta,-2/sigsen^2,2/sigsen)-normcdf(0,-2/sigsen^2,2/sigsen));
            x(itrl,inlo) = x(itrl,inlo)+ ...
                log(normcdf(0,+2/sigsen^2,2/sigsen)-normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(0,-2/sigsen^2,2/sigsen)-normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add inference noise => ALWAYS
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
        end
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_selbias(h,siginf,sigsel,beta)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_selepsi(h,siginf,sigsel,epsi)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

    function [c] = simc_idindep(h,siginf,sigsel,delta1,delta2,delta3)
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
        % compute reversal metrics
        c = getc(dat.t,dat.c,r,dat.w);
    end

end

function [y] = update_prior(x,h)
% update prior belief wrt hazard rate
y = x+log((1-h)./h+exp(-x))-log((1-h)./h+exp(+x));
end
