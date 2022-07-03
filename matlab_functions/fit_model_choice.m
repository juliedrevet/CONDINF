function [out] = fit_model_choice(cfg)
%  FIT_MODEL_CHOICE  Fit model to choice data
%
%  This function fits a desired model to reversal and repetition curves:
%    * modtype = 'flatinf' => flat inference
%    * modtype = 'hrchinf' => hierarchical inference
%    * modtype = 'senbias' => perceptual bias at sensory stage
%    * modtype = 'inflaps' => inference lapses
%    * modtype = 'infdisc' => conditional inference
%
%  Set faltype = 1 for explicit false alarms (visible in dat.s)
%  or faltype = 2 for implicit false alarms (erased from dat.s and inferred
%  online by the particle filter).
%
%  This function requires VBMC v1.0 (https://github.com/lacerbi/vbmc/).
%
%  Valentin Wyart <valentin.wyart@ens.fr> - June 2020
%
%  NB: Julie Drevet: The selection bias and selection lapses are parameters by default.
%  model added as modtype to ease function use:
%    * modtype = 'selbias' => repetition bias at selection stage

% check presence of required parameters
if ~all(isfield(cfg,{'t','s','r','pfal'}))
    error('Incomplete data information!');
end
if ~isfield(cfg,'modtype')
    error('Missing model type!');
end
% check presence of optional parameters
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
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

% create data structure
dat   = [];
dat.t = cfg.t; % trial number in block
dat.s = cfg.s; % stimulus (+1/-1)
dat.r = cfg.r; % response (+1/-1)

% get model type
modtype = lower(cfg.modtype);
if ~ismember(modtype,{'flatinf','hrchinf','senbias','selbias','inflaps','infdisc','infwght'})
    error('Undefined model type!');
end

if strcmp(modtype,'hrchinf')
    modtype = 'senbias';
    cfg.lambda = 0;
end

% selection bias parameter per default in all models, equivalent to fit
% perceptual bias model without lambda
if strcmp(modtype,'selbias') 
    modtype = 'senbias';
    cfg.lambda = 0;
end

ntrl = numel(dat.t); % number of trials
nsmp = cfg.nsmp; % number of samples
nres = cfg.nres; % number of re-samples for bootstrapping
nval = cfg.nval; % number of function evaluations for validation
verbose = cfg.verbose; % display level

% define response observations
isrpos = dat.r > 0;
isrneg = dat.r < 0;

% compute sensory noise
pfal = cfg.pfal; % false-alarm rate
% set sensory noise
sigsen = max(1/norminv(1-pfal),1e-6);


% define model parameters
pnam = {}; % name
pmin = []; % minimum value
pmax = []; % maximum value
pfun = {}; % log-prior function
pini = []; % initial value
pplb = []; % plausible lower bound
ppub = []; % plausible upper bound
% perceived hazard rate ~ beta(2,18)
pnam{1,1} = 'h';
pmin(1,1) = 0;
pmax(1,1) = 0.5;
pfun{1,1} = @(x)betapdf(x,2,18);
pini(1,1) = betastat(2,18);
pplb(1,1) = betainv(0.1587,2,18);
ppub(1,1) = betainv(0.8413,2,18);
% inference noise ~ gamma(4,1/4)
pnam{1,2} = 'siginf';
pmin(1,2) = 1e-6;
pmax(1,2) = 5;
pfun{1,2} = @(x)gampdf(x,4,0.25);
pini(1,2) = gamstat(4,0.25);
pplb(1,2) = gaminv(0.1587,4,0.25);
ppub(1,2) = gaminv(0.8413,4,0.25);
% selection noise ~ gamma(4,1/4)
pnam{1,3} = 'sigsel';
pmin(1,3) = 1e-6;
pmax(1,3) = 5;
pfun{1,3} = @(x)gampdf(x,4,0.25);
pini(1,3) = gamstat(4,0.25);
pplb(1,3) = gaminv(0.1587,4,0.25);
ppub(1,3) = gaminv(0.8413,4,0.25); 
% selection lapses ~ beta(1,39) 
pnam{1,4} = 'epsi';
pmin(1,4) = 1e-3; 
pmax(1,4) = 1;
pfun{1,4} = @(x)betapdf(x,1,39);
pini(1,4) = betastat(1,39); 
pplb(1,4) = betainv(0.1587,1,39); 
ppub(1,4) = betainv(0.8413,1,39); 
% selection bias ~ gamma(1,1/4)
pnam{1,5} = 'beta';
pmin(1,5) = 0;
pmax(1,5) = 5; 
pfun{1,5} = @(x)gampdf(x,1,0.25);
pini(1,5) = gamstat(1,0.25);
pplb(1,5) = gaminv(0.1587,1,0.25); 
ppub(1,5) = gaminv(0.8413,1,0.25); 
switch modtype
    case 'senbias'
        % sensory bias strength ~ gamma(1,1)
        pnam{1,6} = 'lambda';
        pmin(1,6) = 0;
        pmax(1,6) = 5;
        pfun{1,6} = @(x)gampdf(x,1,1);
        pini(1,6) = gamstat(1,1);
        pplb(1,6) = gaminv(0.1587,1,1);
        ppub(1,6) = gaminv(0.8413,1,1);
    case 'inflaps'
        % fraction of inference lapses ~ beta(1,3) 
        pnam{1,6} = 'plaps';
        pmin(1,6) = 0;
        pmax(1,6) = 1;
        pfun{1,6} = @(x)betapdf(x,1,3);
        pini(1,6) = betastat(1,3);
        pplb(1,6) = betainv(0.1587,1,3);
        ppub(1,6) = betainv(0.8413,1,3);
    case {'infdisc','infwght'}
        % reliability threshold ~ gamma(1,1)
        pnam{1,6} = 'delta';
        pmin(1,6) = 0;
        pmax(1,6) = 5;
        pfun{1,6} = @(x)gampdf(x,1,1);
        pini(1,6) = gamstat(1,1);
        pplb(1,6) = gaminv(0.1587,1,1);
        ppub(1,6) = gaminv(0.8413,1,1);
end

% set number of parameters
npar = numel(pnam);

% define fixed parameters
pfix = cell(1,npar);
for i = 1:npar
    if isfield(cfg,pnam{i})
        % clip parameter value inside valid range
        cfg.(pnam{i}) = min(max(cfg.(pnam{i}),pmin(i)),pmax(i));
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
    out.modtype = cfg.modtype;
    
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

% store filtered response probabilities
[out.p,out.x] = getp(phat{:});

% store configuration structure
out.cfg = cfg;
    
if nfit > 0 && nval > 0
    % perform function evaluations for validation
    l_val = nan(nval,1); % estimates for each validation sample
    l_sd_val = nan(nval,1); % bootstrap s.d. estimates for each validation sample
    for ival = 1:nval
        [l_val(ival),l_sd_val(ival)] = fun(xhat);
    end
    out.l_val = l_val; % estimates
    out.l_sd_val = l_sd_val; % bootstrap s.d. estimates
end

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
        % compute filtered response probability
        p = getp(varargin{:});
        % compute bootstrapped log-likelihood s.d.
        lres = nan(nres,1);
        for ires = 1:nres
            jres = randsample(nsmp,nsmp,true);
            pres = mean(p(:,jres),2);
            lres(ires) = ...
                sum(log(pres(isrpos)))+ ...
                sum(log(1-pres(isrneg)));
        end
        ll_sd = max(std(lres),1e-6);
        % compute log-likelihood
        p = mean(p,2);
        ll = ...
            sum(log(p(isrpos)))+ ...
            sum(log(1-p(isrneg)));
    end

    function [p,x] = getp(varargin)
        switch modtype
            case 'flatinf'
                [p,x] = getp_flatinf(varargin{:});
            case 'senbias'
                [p,x] = getp_senbias(varargin{:});
            case 'inflaps'
                [p,x] = getp_inflaps(varargin{:});
            case 'infdisc'
                [p,x] = getp_infdisc(varargin{:});
            case 'infwght'
                [p,x] = getp_infwght(varargin{:});
        end
    end

    function [p,x] = getp_flatinf(h,siginf,sigsel,epsi,beta)
        % merge sensory noise into inference noise
        siginf = sqrt(4/sigsen^2+siginf^2);
        % run particle filter
        p = zeros(ntrl,nsmp); % filtered response probability
        x = zeros(ntrl,nsmp); % filtered log-belief
        for itrl = 1:ntrl
            % compute response probability
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
                x_c = 0;
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
                x_c = -beta*dat.r(itrl-1);
            end
            x_t = x(itrl,:)+dat.s(itrl)*2/sigsen^2;
            p(itrl,:) = 1-normcdf(x_c,x_t,sqrt(siginf^2+sigsel^2));
            p(itrl,:) = (1-epsi)*p(itrl,:)+epsi*0.5;
            % flag epsilon trials
            ie = rand(1,nsmp) < epsi;
            if any(ie)
                % update log-belief for epsilon trials
                xs = x_t+siginf*randn(1,nsmp);
                x(itrl,ie) = xs(ie);
            end
            % filter log-belief
            x_f = tnormrnd(x_t,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            x_t = normrnd( ...
                (sigsel^2*x_t+siginf^2*x_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            % update log-belief for regular trials
            xs = x_t;
            x(itrl,~ie) = xs(~ie);
        end
    end

    function [p,x] = getp_senbias(h,siginf,sigsel,epsi,beta,lambda)
        % run particle filter
        p = zeros(ntrl,nsmp); % filtered response probability
        x = zeros(ntrl,nsmp); % filtered log-belief
        for itrl = 1:ntrl
            % compute response probability
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
                x_c = 0;
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
                x_c = -beta*dat.r(itrl-1);
            end
            ph0 = 1-normcdf(-x(itrl,:)*dat.s(itrl)*lambda,2/sigsen^2,2/sigsen);
            pf0 = 1-ph0;
            xp = ...
                log(1-normcdf(x(itrl,:)*lambda,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(x(itrl,:)*lambda,-2/sigsen^2,2/sigsen));
            xn = ...
                log(normcdf(x(itrl,:)*lambda,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(x(itrl,:)*lambda,-2/sigsen^2,2/sigsen));
            if dat.s(itrl) > 0
                xh = x(itrl,:)+xp;
                xf = x(itrl,:)+xn;
            else
                xh = x(itrl,:)+xn;
                xf = x(itrl,:)+xp;
            end
            ph = 1-normcdf(x_c,xh,sqrt(siginf^2+sigsel^2));
            pf = 1-normcdf(x_c,xf,sqrt(siginf^2+sigsel^2));
            p(itrl,:) = ph.*ph0+pf.*pf0;
            p(itrl,:) = (1-epsi)*p(itrl,:)+epsi*0.5;
            % flag epsilon trials
            ie = rand(1,nsmp) < epsi;
            if any(ie)
                % update log-belief for epsilon trials
                [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph0;pf0]./(ph0+pf0),1)));
                xs = [xh;xf];
                xs = xs+siginf*randn(2,nsmp);
                xs = xs(sub2ind([2,nsmp],is,1:nsmp));
                x(itrl,ie) = xs(ie);
            end
            % filter log-belief
            xh_f = tnormrnd(xh,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xf_f = tnormrnd(xf,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xh = normrnd( ...
                (sigsel^2*xh+siginf^2*xh_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            xf = normrnd( ...
                (sigsel^2*xf+siginf^2*xf_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            if dat.r(itrl) < 0
                ph = 1-ph;
            end
            ph = ph.*ph0;
            if dat.r(itrl) < 0
                pf = 1-pf;
            end
            pf = pf.*pf0;
            % update log-belief for regular trials
            [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph;pf]./(ph+pf),1)));
            xs = [xh;xf];
            xs = xs(sub2ind([2,nsmp],is,1:nsmp));
            x(itrl,~ie) = xs(~ie);
        end
    end

    function [p,x] = getp_inflaps(h,siginf,sigsel,epsi,beta,plaps)
        % compute base rates
        ph0 = (1-plaps)*(1-pfal); % base hit rate
        pf0 = (1-plaps)*pfal; % base false-alarm rate
        pl0 = plaps; % base lapse rate
        % run particle filter
        p = zeros(ntrl,nsmp); % filtered response probability
        x = zeros(ntrl,nsmp); % filtered log-belief
        for itrl = 1:ntrl
            % compute response probability
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
                x_c = 0;
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
                x_c = -beta*dat.r(itrl-1);
            end
            xh = x(itrl,:)+dat.s(itrl)*log((1-pfal)/pfal);
            xf = x(itrl,:)-dat.s(itrl)*log((1-pfal)/pfal);
            xl = x(itrl,:);
            ph = 1-normcdf(x_c,xh,sqrt(siginf^2+sigsel^2));
            pf = 1-normcdf(x_c,xf,sqrt(siginf^2+sigsel^2));
            pl = xl > 0;
            p(itrl,:) = ph.*ph0+pf.*pf0+pl.*pl0;
            p(itrl,:) = (1-epsi)*p(itrl,:)+epsi*0.5;
            % flag epsilon trials
            ie = rand(1,nsmp) < epsi;
            if any(ie)
                % update log-belief for epsilon trials
                [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph0;pf0;pl0]./(ph0+pf0+pl0),1)));
                xs = [xh;xf;xl];
                xs([1,2],:) = xs([1,2],:)+siginf*randn(2,nsmp);
                xs = xs(sub2ind([3,nsmp],is,1:nsmp));
                x(itrl,ie) = xs(ie);
            end
            % filter log-belief
            xh_f = tnormrnd(xh,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xf_f = tnormrnd(xf,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xh = normrnd( ...
                (sigsel^2*xh+siginf^2*xh_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            xf = normrnd( ...
                (sigsel^2*xf+siginf^2*xf_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            if dat.r(itrl) < 0
                ph = 1-ph;
            end
            ph = ph*ph0;
            if dat.r(itrl) < 0
                pf = 1-pf;
            end
            pf = pf*pf0;
            if dat.r(itrl) < 0
                pl = 1-pl;
            end
            pl = pl*pl0;
            % update log-belief for regular trials
            [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph;pf;pl]./(ph+pf+pl),1)));
            xs = [xh;xf;xl];
            xs = xs(sub2ind([3,nsmp],is,1:nsmp));
            x(itrl,~ie) = xs(~ie);
        end
    end

    function [p,x] = getp_infdisc(h,siginf,sigsel,epsi,beta,delta)
        % compute base rates
        ph0 = 1-normcdf(+delta,2/sigsen^2,2/sigsen); % base hit rate
        pf0 = normcdf(-delta,2/sigsen^2,2/sigsen); % base false-alarm rate
        pl0 = 1-ph0-pf0; % base lapse rate
        % run particle filter
        p = zeros(ntrl,nsmp); % filtered response probability
        x = zeros(ntrl,nsmp); % filtered log-belief
        for itrl = 1:ntrl
            % compute response probability
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
                x_c = 0;
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
                x_c = -beta*dat.r(itrl-1);
            end
            xh = x(itrl,:)+dat.s(itrl)*( ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen)));
            xf = x(itrl,:)-dat.s(itrl)*( ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen)));
            xl = x(itrl,:);
            ph = 1-normcdf(x_c,xh,sqrt(siginf^2+sigsel^2));
            pf = 1-normcdf(x_c,xf,sqrt(siginf^2+sigsel^2));
            pl = xl > 0;
            p(itrl,:) = ph.*ph0+pf.*pf0+pl.*pl0;
            p(itrl,:) = (1-epsi)*p(itrl,:)+epsi*0.5;
            % flag epsilon trials
            ie = rand(1,nsmp) < epsi;
            if any(ie)
                % update log-belief for epsilon trials
                [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph0;pf0;pl0]./(ph0+pf0+pl0),1)));
                xs = [xh;xf;xl];
                xs([1,2],:) = xs([1,2],:)+siginf*randn(2,nsmp);
                xs = xs(sub2ind([3,nsmp],is,1:nsmp));
                x(itrl,ie) = xs(ie);
            end
            % filter log-belief
            xh_f = tnormrnd(xh,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xf_f = tnormrnd(xf,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xh = normrnd( ...
                (sigsel^2*xh+siginf^2*xh_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            xf = normrnd( ...
                (sigsel^2*xf+siginf^2*xf_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            if dat.r(itrl) < 0
                ph = 1-ph;
            end
            ph = ph*ph0;
            if dat.r(itrl) < 0
                pf = 1-pf;
            end
            pf = pf*pf0;
            if dat.r(itrl) < 0
                pl = 1-pl;
            end
            pl = pl*pl0;
            % update log-belief for regular trials
            [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph;pf;pl]./(ph+pf+pl),1)));
            xs = [xh;xf;xl];
            xs = xs(sub2ind([3,nsmp],is,1:nsmp));
            x(itrl,~ie) = xs(~ie);
        end
    end

    function [p,x] = getp_infwght(h,siginf,sigsel,epsi,beta,delta)
        % compute base rates
        ph0 = 1-normcdf(+delta,2/sigsen^2,2/sigsen); % base hit rate
        pf0 = normcdf(-delta,2/sigsen^2,2/sigsen); % base false-alarm rate
        pl0 = 1-ph0-pf0; % base lapse rate
        % run particle filter
        p = zeros(ntrl,nsmp); % filtered response probability
        x = zeros(ntrl,nsmp); % filtered log-belief
        for itrl = 1:ntrl
            % compute response probability
            if dat.t(itrl) == 1
                x(itrl,:) = 1e-6*sign(randn(1,nsmp));
                x_c = 0;
            else
                x(itrl,:) = update_prior(x(itrl-1,:),h);
                x_c = -beta*dat.r(itrl-1);
            end
            xh = x(itrl,:)+dat.s(itrl)*( ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen)));
            xf = x(itrl,:)-dat.s(itrl)*( ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen)));
            xl = x(itrl,:);
            ph = 1-normcdf(x_c,xh,sqrt(siginf^2+sigsel^2));
            pf = 1-normcdf(x_c,xf,sqrt(siginf^2+sigsel^2));
            pl = xl > 0;
            p(itrl,:) = ph.*ph0+pf.*pf0+pl.*pl0;
            p(itrl,:) = (1-epsi)*p(itrl,:)+epsi*0.5;
            % flag epsilon trials
            ie = rand(1,nsmp) < epsi;
            if any(ie)
                % update log-belief for epsilon trials
                [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph0;pf0;pl0]./(ph0+pf0+pl0),1)));
                xs = [xh;xf;xl];
                xs([1,2],:) = xs([1,2],:)+siginf*randn(2,nsmp);
                xs = xs(sub2ind([3,nsmp],is,1:nsmp));
                x(itrl,ie) = xs(ie);
            end
            % filter log-belief
            xh_f = tnormrnd(xh,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xf_f = tnormrnd(xf,sqrt(siginf^2+sigsel^2),x_c,dat.r(itrl));
            xh = normrnd( ...
                (sigsel^2*xh+siginf^2*xh_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            xf = normrnd( ...
                (sigsel^2*xf+siginf^2*xf_f)./(siginf^2+sigsel^2), ...
                sqrt(siginf^2*sigsel^2/(siginf^2+sigsel^2)));
            if dat.r(itrl) < 0
                ph = 1-ph;
            end
            ph = ph*ph0;
            if dat.r(itrl) < 0
                pf = 1-pf;
            end
            pf = pf*pf0;
            if dat.r(itrl) < 0
                pl = 1-pl;
            end
            pl = pl*pl0;
            % update log-belief for regular trials
            [~,is] = max(bsxfun(@lt,rand(1,nsmp),cumsum([ph;pf;pl]./(ph+pf+pl),1)));
            xs = [xh;xf;xl];
            xs = xs(sub2ind([3,nsmp],is,1:nsmp));
            x(itrl,~ie) = xs(~ie);
        end
    end

end

function [x] = tnormrnd(m,s,t,d)
% sample from truncated normal distribution
if d == 1 % positive truncation
    x = +rpnormv(+m-t,s)+t;
else % negative truncation
    x = -rpnormv(-m+t,s)+t;
end
end

function [y] = update_prior(x,h)
% update prior belief as a function of hazard rate
y = x+log((1-h)./h+exp(-x))-log((1-h)./h+exp(+x));
end
