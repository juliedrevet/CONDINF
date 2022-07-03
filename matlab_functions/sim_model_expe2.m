function [out] = sim_model_expe2(cfg)
%  SIM_MODEL  Simulate model
%
%  This function simulates a desired model:
%    * modtype = 'flatinf' => flat inference
%    * modtype = 'nonstab' => hierarchical inference w/o stabilization
%    * modtype = 'senbias' => sensory bias
%    * modtype = 'inflaps' => inference lapses
%    * modtype = 'infdisc' => inference discard
%    * modtype = 'selbias' => selection bias
%    * modtype = 'selepsi' => epsilon-greedy selection
%    * modtype = 'infwght' => weighted inference (added after review)
%
%  Valentin Wyart <valentin.wyart@ens.fr> - June 2020/ Julie Drevet 2021

if ~isfield(cfg,'seed') % seed for random number generator
    cfg.seed = round(sum(100*clock)); % seed based on current time
end
% initialize random number generator
sd = RandStream('mt19937ar','Seed',cfg.seed);
RandStream.setGlobalStream(sd);

% define useful function handles
logit = @(x)log(x./(1-x)); % logit function

% create data structure
dat   = [];
dat.t = cfg.t; % trial number in block
dat.s = cfg.s; % stimulus: will be normalized according to strength
dat.c = cfg.s; % stimulus: as underlying category (+1/-1)
dat.w = cfg.w; % stimulus weight (1/2/3)

ntrl = numel(dat.t); % number of trials
nsmp = cfg.nsmp; % number of simulations

% compute sensory noise
pfal   = cfg.pfal; % false-alarm rate
sigsen = 1/norminv(1-pfal(2)); % sensory noise
strsen = norminv(1-pfal)/norminv(1-pfal(2)); % sensory strengths
% normalize stimulus
for iw = 1:3
    dat.s(dat.w == iw) = strsen(iw)*dat.s(dat.w == iw);
end

% get model type
modtype = lower(cfg.modtype);
if ~ismember(modtype,{'flatinf','nonstab','senbias', ...
        'inflaps','infdisc','selbias','selepsi','infwght'})
    error('Undefined model type!');
end
if strcmp(modtype,'nonstab')
    modtype = 'infdisc';
    cfg.delta = 0;
end

cfg.h = min(max(cfg.h,1e-9),0.5);
% get model variables
var = {cfg.h,cfg.siginf,cfg.sigsel};
switch modtype
    case 'senbias', var{end+1} = cfg.lambda;
    case 'inflaps', var{end+1} = cfg.plaps;
    case 'infdisc', var{end+1} = cfg.delta;
    case 'selbias', var{end+1} = cfg.beta;
    case 'selepsi', var{end+1} = cfg.epsi;
    case 'infwght', var{end+1} = cfg.delta;
end

% simulate responses
r = simr(var{:});

% compute reversal metrics
c = getc(dat.t,dat.c,r,dat.w);

% get response repetitions
rep = [nan(1,nsmp);r(2:end,:) == r(1:end-1,:)];
rep(dat.t == 1,:) = nan;

% create output structure
out = [];
out.r = r; % simulated responses
out.c = c; % reversal metrics
out.pcor = mean(bsxfun(@eq,r,dat.c),[1,2]); % response accuracy
out.prep = mean(nanmean(rep,1),2); % response stability

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
            case 'infwght' % weighted inference
                r = simr_infwght(varargin{:});
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
            r(itrl,r(itrl,:)==0) = sign(randn(1,nnz(r(itrl,:)==0)));
        end
    end

    function [r] = simr_infwght(h,siginf,sigsel,delta)
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
            iphi = s > +delta;
            inhi = s < -delta;
            iplo = s < +delta & s > 0; % below threshold
            inlo = s > -delta & s < 0;
            % add log-evidence
            x(itrl,iphi) = x(itrl,iphi)+ ...
                log(1-normcdf(+delta,+2/sigsen^2,2/sigsen))- ...
                log(1-normcdf(+delta,-2/sigsen^2,2/sigsen));
            x(itrl,inhi) = x(itrl,inhi)+ ...
                log(normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add log-evidence BELOW THRESHOLD
            x(itrl,iplo) = x(itrl,iplo)+ ...
                log(normcdf(+delta,+2/sigsen^2,2/sigsen)-normcdf(0,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(+delta,-2/sigsen^2,2/sigsen)-normcdf(0,-2/sigsen^2,2/sigsen));
            x(itrl,inlo) = x(itrl,inlo)+ ...
                log(normcdf(0,+2/sigsen^2,2/sigsen)-normcdf(-delta,+2/sigsen^2,2/sigsen))- ...
                log(normcdf(0,-2/sigsen^2,2/sigsen)-normcdf(-delta,-2/sigsen^2,2/sigsen));
            % add inference noise ALWAYS
            x(itrl,:) = x(itrl,:)+siginf*randn(1,nsmp);
            % sample response with selection noise
            r(itrl,:) = sign(x(itrl,:)+sigsel*randn(1,nsmp));
            r(itrl,r(itrl,:)==0) = sign(randn(1,nnz(r(itrl,:)==0)));
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

end

function [y] = update_prior(x,h)
% update prior belief wrt hazard rate
y = x+log((1-h)./h+exp(-x))-log((1-h)./h+exp(+x));
end
