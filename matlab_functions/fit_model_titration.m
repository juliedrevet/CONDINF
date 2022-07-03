function [out] = fit_model_titration(cfg)
%  FIT_MODEL_TITR  Fit confidence model to titration data
%
%  Usage: [out] = fit_model_titr(cfg)
%
%  The model is controlled by three free parameters:
%    * sigsen = sensory noise
%    * sigcon = confidence noise
%    * thrcon = confidence threshold
%
%  Valentin Wyart <valentin.wyart@ens.fr>

% check configuration structure
if ~all(isfield(cfg,{'stim','resp','conf'}))
    error('Missing experiment data!');
end
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
end
if ~isfield(cfg,'nres')
    cfg.nres = 1e3;
end
if ~isfield(cfg,'verbose')
    cfg.verbose = false;
end

stim = cfg.stim; % stimulus (+1:light/-1:dark)
resp = cfg.resp; % response (+1:light/-1:dark)
conf = cfg.conf; % confidence (+1:high/-1:low)

nsmp = cfg.nsmp; % number of samples
nres = cfg.nres; % number of bootstrap resamples
ntrl = numel(stim); % number of trials

peps = 1e-6; % minimum filtered probability
verbose = cfg.verbose; % verbose VBMC output?

% get factorized response/confidence index
ifact = sub2ind([ntrl,4],(1:ntrl)',2*(resp < 0)+(conf < 0)+1);

% define model parameters
pnam = {}; % name
pmin = []; % minimum value
pmax = []; % maximum value
pfun = {}; % log-prior function
pini = []; % initial value
pplb = []; % plausible lower bound
ppub = []; % plausible upper bound
% 1/ sensory noise ~ gamma(a,b)
m = 1/norminv(0.8); % mean
a = sqrt(m)*32; % shape
b = sqrt(m)/32; % scale
pnam{1,1} = 'sigsen';
pmin(1,1) = 0.001;
pmax(1,1) = 5;
pfun{1,1} = @(x)gampdf(x,a,b);
pini(1,1) = gamstat(a,b);
pplb(1,1) = gaminv(0.15,a,b);
ppub(1,1) = gaminv(0.85,a,b);
% 2/ confidence noise ~ gamma(a,b)
a = 1; % shape
b = 2; % scale
pnam{1,2} = 'sigcon';
pmin(1,2) = 0.001;
pmax(1,2) = 10;
pfun{1,2} = @(x)gampdf(x,a,b);
pini(1,2) = gamstat(a,b);
pplb(1,2) = gaminv(0.15,a,b);
ppub(1,2) = gaminv(0.85,a,b);
% 3/ confidence threshold ~ normal(m,s)
m = 1; % mean
s = 2; % standard deviation
pnam{1,3} = 'thrcon';
pmin(1,3) = -10;
pmax(1,3) = +10;
pfun{1,3} = @(x)normpdf(x,m,s);
pini(1,3) = normstat(m,s);
pplb(1,3) = norminv(0.15,m,s);
ppub(1,3) = norminv(0.85,m,s);

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
        options.Display = 'none';
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

% store configuration structure
out.cfg = cfg;

% store filtered response/confidence probabilities
[out.presp,out.pconf] = getp(phat{:});
out.pconf = mean(out.pconf,2); % average across samples

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
        % get filtered response/confidence probabilities
        [presp,pconf] = getp(varargin{:});
        % compute bootstrap estimate of log-likelihood s.d.
        ll_res = nan(nres,1);
        for ires = 1:nres
            jres = randsample(nsmp,nsmp,true);
            pconf_res = mean(pconf(:,jres),2);
            pfact_res = cat(2, ...
                presp.*pconf_res,presp.*(1-pconf_res), ...
                (1-presp).*pconf_res,(1-presp).*(1-pconf_res));
            ll_res(ires) = sum(log(max(pfact_res(ifact),peps)));
        end
        ll_sd = max(std(ll_res),1e-6);
        % compute log-likelihood estimate
        pconf = mean(pconf,2);
        pfact = cat(2, ...
            presp.*pconf,presp.*(1-pconf), ...
            (1-presp).*pconf,(1-presp).*(1-pconf));
        ll = sum(log(max(pfact(ifact),peps)));
    end

    function [presp,pconf] = getp(sigsen,sigcon,thrcon)
        % compute filtered response/confidence probabilities
        presp = nan(ntrl,1); % response probabilities
        pconf = nan(ntrl,nsmp); % confidence probabilities
        for itrl = 1:ntrl
            presp(itrl) = 1-normcdf(0,stim(itrl),sigsen);
            xt = tnormrnd(stim(itrl(ones(1,nsmp))),sigsen,0,resp(itrl));
            pconf(itrl,:) = 1-normcdf(thrcon,abs(xt),sigcon);
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