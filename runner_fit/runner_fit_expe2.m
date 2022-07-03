%% expe2_runner_fit
% needed on the path for fitting models to participants' behaviour:
%       *VBMC v1.0 (https://github.com/lacerbi/vbmc/)
%       *BADS v1.0 (https://github.com/lacerbi/bads/)

% clear workspace
clearvars
close all
clc

addpath '../matlab_functions';

datapath = '../data_table/';
savepath = '../results/';
if ~exist(savepath,'dir')
    mkdir(savepath);
end
%% check participants 
nstd = 3.5;
pcor = [];
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');
% get list of participants
participants = unique(dat.participant);
ndata = numel(participants);
% compute performance
for idata = 1:ndata
    ifilt = dat.participant == participants(idata);
    s = dat.stim(ifilt);
    r = dat.resp(ifilt);
    pcor = cat(1,pcor,nanmean(r == s));
end
% compute mean+std performance leaving-one-out participant
pcor_loo_avg = nan(ndata,1);
pcor_loo_std = nan(ndata,1);
for idata = 1:ndata
    pcor_loo_avg(idata) = mean(pcor(setdiff(1:ndata,idata)));
    pcor_loo_std(idata) = std(pcor(setdiff(1:ndata,idata)));
end
I = find((pcor<(pcor_loo_avg+pcor_loo_std*(-nstd))) | (pcor > (pcor_loo_avg+pcor_loo_std*(+nstd))));
if ~isempty(I)
    fprintf('remove participant S%02d from analysis\n',participants(I))
end
participants = setdiff(unique(dat.participant),I);
ndata = numel(participants);

%% FIT NOISY INFERENCE MODEL
fprintf('<strong>Fitting noisy inference model</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);       % trial number in block
    s = dat.stim(ifilt);        % stimulus (+1/-1)
    w = dat.wght(ifilt);        % stimulus weight or strength (1/2/3)
    r = dat.resp(ifilt);        % response (+1/-1)

    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;               % trial number in block
    cfg_fit.s        = s;               % stimulus
    cfg_fit.w        = w;               % strength
    cfg_fit.r        = r;               % response
    cfg_fit.nsmp     = 1e3;             % number of samples
    cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
    cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
    cfg_fit.sigsel   = 0;               % no selection noise
    cfg_fit.modtype  = 'nonstab';
    
    out_fit{idata} = fit_model_revrepwgh(cfg_fit);
end
% save results
save(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit');

%% FIT STABILIZING MODELS
fprintf('<strong>Fitting stabilizing models</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

models = {'senbias','inflaps','infdisc','selbias','selepsi'};
nmodel = numel(models);
out_fit = cell(ndata,nmodel);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    w = dat.wght(ifilt);    % weight or strength (1/2/3)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;               % trial number in block
    cfg_fit.s        = s;               % stimulus
    cfg_fit.w        = w;               % strength
    cfg_fit.r        = r;               % response
    cfg_fit.nsmp     = 1e3;             % number of samples
    cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
    cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
    cfg_fit.sigsel   = 0;               % no selection noise
    
    for imodel = 1:nmodel
        % fit specified model
        fprintf('\tmodel %s...',models{imodel});
        cfg_tmp = cfg_fit;
        cfg_tmp.modtype = models{imodel};
        out_fit{idata,imodel} = fit_model_revrepwgh(cfg_tmp);
    end
end
% save results
save(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');

%% FIT ONE RELIABILITY THRESHOLD PER STIMULUS STRENGTH
fprintf('<strong>Fitting one reliability threshold per stimulus strength for conditional inference model</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    w = dat.wght(ifilt);    % strength (1/2/3)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;               % trial number in block
    cfg_fit.s        = s;               % stimulus
    cfg_fit.w        = w;               % strength
    cfg_fit.r        = r;               % response
    cfg_fit.nsmp     = 1e3;             % number of samples
    cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
    cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
    cfg_fit.sigsel   = 0;               % no selection noise
    cfg_fit.modtype  = 'idindep';
    
    out_fit{idata} = fit_model_revrepwgh(cfg_fit);
end
% save results
save(fullfile(savepath,'expe2_results_fit_conditional_inference_independent_thresholds.mat'),'out_fit');

%% FIT CONFIDENCE MODEL TO TITRATION DATA
fprintf('<strong>Fitting confidence model to titration data</strong>\n\n')

clearvars -except datapath savepath ndata participants
% load titration results
load(fullfile(datapath,'expe2_titration.mat'));
dat = expe2_titration;
clear('expe2_titration');

% set analysis parameters
nburnin = 50; % number of burnin trials at the beginning of each session

% get list of participants
participants = unique(dat.participant);
ndata = numel(participants);

% initialize output
out_fit = cell(ndata,1); % fitting output

for idata = 1:ndata
    fprintf('\tParticipant %02d/%02d...',idata,ndata)
    % get data
    sess = []; % session number (1/2)
    strc = []; % staircase side (+1:light/-1:dark)
    prop = []; % stimulus proportion of light
    resp = []; % response (+1:light/-1:dark)
    conf = []; % confidence (+1:high/-1:low)
    % define conditions
    icond = { ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block > 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block > 8)};
    ncond = numel(icond);
    for i = 1:ncond
        % get condition
        ifilt = icond{i};
        % exclude burnin trials
        ifilt = ifilt(nburnin+1:end);
        nfilt = numel(ifilt);
        % get data
        sess = cat(1,sess,ceil(dat.block(ifilt)/8));
        strc = cat(1,strc,dat.staircase(ifilt));
        prop = cat(1,prop,dat.prop(ifilt));
        resp = cat(1,resp,dat.resp(ifilt));
        conf = cat(1,conf,dat.conf(ifilt));
    end
    
    % fit confidence model
    cfg      = [];
    cfg.stim = strc;
    cfg.resp = resp;
    cfg.conf = conf;
    cfg.nsmp = 1e2; % number of particle samples
    cfg.nres = 1e2; % number of bootstrap resamples
    out = fit_model_titration(cfg);
    out_fit{idata} = out;
    fprintf('done\n');
end
% save output
save(fullfile(savepath,'expe2_results_fit_confidence.mat'),'out_fit');

%% SUPPLEMENTARY: FIT STABILIZING PARAMETERS BEST_FITTING PARTICIPANTS' OVERALL SWITCH RATE
fprintf('<strong>Fitting stabilizing parameters best-fitting participants'' overall switch rate</strong>\n')

clearvars -except datapath savepath ndata participants

% start from unstabilized best-fitting model parameters and add the
% stabilizing parameter
load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit')
out_fit_noisy = out_fit;
models = {'senbias','inflaps','infdisc','selbias','selepsi'};
nmodel = numel(models);
out_fit = cell(ndata,nmodel);

for idata = 1:ndata
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    cfg_fit         = out_fit_noisy{idata}.cfg;
    cfg_fit.h       = out_fit_noisy{idata}.h;
    cfg_fit.siginf  = out_fit_noisy{idata}.siginf;
    cfg_fit.sigsel  = 0;
    cfg_fit.ftype   = 'negllh';
    cfg_fit.verbose = false;
    for imodel = 1:nmodel
        fprintf('\tmodel %s...',models{imodel});
        cfg_fit.modtype = models{imodel}; 
        out_fit{idata,imodel} = fit_overall_stab_expe2(cfg_fit);
        fprintf('done.\n')
    end
end
% save results
save(fullfile(savepath,'expe2_results_fit_stabilizing_parameters.mat'),'out_fit');

%% SUPPLEMENTARY: EX-ANTE RECOVERY STABILIZING MODELS
fprintf('<strong>Ex-ante recovery between stabilizing models</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

models = {'senbias','inflaps','infdisc','selbias','selepsi'};
nmodel = numel(models);
nsim   = 10; % number of simulations to be recovered

elbo = nan(ndata,nmodel,nmodel,nsim);

% simulation parameters
% 1/ perceived hazard rate
h0 = betastat(2,18);
% 2/ inference noise
siginf0 = gamstat(4,0.25);
% 3/ selection noise => no selection noise
% 4/ stab param
lambda0 = gamstat(1,0.75);
plaps0  = betastat(1,9);
delta0  = gamstat(1,0.75);
beta0   = gamstat(1,0.25);
epsi0   = betastat(1,9);

for isim = 1:nsim
    for idata = 1:ndata
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);	% trial number in block
        s = dat.stim(ifilt);	% stimulus (+1/-1)
        w = dat.wght(ifilt);    % strength (1/2/3)
        for imodel = 1:nmodel
            fprintf('processing participant %02d/%02d - simulation model %s - simulation %02d\n\n',idata,ndata,models{imodel},isim);
            modtype = models{imodel};
            
            % 1/ simulate model
            cfg_sim = [];
            cfg_sim.modtype = modtype;
            cfg_sim.t = t;                      
            cfg_sim.s = s;                      
            cfg_sim.w = w;
            cfg_sim.pfal = [0.30 0.20 0.10];
            cfg_sim.nsmp = 1;
            cfg_sim.seed = isim * idata * 100;
            % simulation parameters
            cfg_sim.h       = h0;
            cfg_sim.siginf  = siginf0;
            cfg_sim.sigsel  = 0;
            switch cfg_sim.modtype
                case 'senbias', cfg_sim.lambda  = lambda0;
                case 'inflaps', cfg_sim.plaps   = plaps0;
                case 'infdisc', cfg_sim.delta   = delta0;
                case 'selbias', cfg_sim.beta    = beta0;
                case 'selepsi', cfg_sim.epsi    = epsi0;
            end
            out = sim_model_expe2(cfg_sim);
            
            % configure model fitting
            cfg_fit          = [];
            cfg_fit.t        = cfg_sim.t;           % trial number in block
            cfg_fit.s        = cfg_sim.s;           % stimulus
            cfg_fit.w        = cfg_sim.w;           % strength
            cfg_fit.r        = out.r;               % response
            cfg_fit.nsmp     = 1e3;                 % number of samples
            cfg_fit.nres     = 1e3;                 % number of re-samples for bootstrapping
            cfg_fit.pfal     = [0.30 0.20 0.10];    % false-alarm rate
            cfg_fit.verbose  = false;               % verbose VBMC output
            cfg_fit.sigsel   = 0;
            
            for imodrec = 1:nmodel
                fprintf('\trecovery with model %s\n\n',models{imodrec});
                cfg_tmp = cfg_fit;
                cfg_tmp.modtype = models{imodrec};
                out_fit = fit_model_revrepwgh(cfg_tmp);
                elbo(idata,imodel,imodrec,isim) = out_fit.elbo;
            end
        end
    end
end
% save results
save(fullfile(savepath,'expe2_results_recovery_stabilizing_models_elbo.mat'),'elbo');

%% SUPPLEMENTARY: FIT CONFIDENCE MODEL TO TITRATION DATA - NO CONFIDENCE NOISE
fprintf('<strong>Fitting confidence model to titration data - no confidence noise</strong>\n\n')

clearvars -except datapath savepath ndata participants
% load titration results
load(fullfile(datapath,'expe2_titration.mat'));
dat = expe2_titration;
clear('expe2_titration');

% set analysis parameters
nburnin = 50; % number of burnin trials at the beginning of each session

% get list of participants
participants = unique(dat.participant);
ndata = numel(participants);

% initialize output
out_fit = cell(ndata,1); % fitting output

for idata = 1:ndata
    fprintf('\tParticipant %02d/%02d...',idata,ndata)
    % get data
    sess = []; % session number (1/2)
    strc = []; % staircase side (+1:light/-1:dark)
    prop = []; % stimulus proportion of light
    resp = []; % response (+1:light/-1:dark)
    conf = []; % confidence (+1:high/-1:low)
    % define conditions
    icond = { ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block > 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block > 8)};
    ncond = numel(icond);
    for i = 1:ncond
        % get condition
        ifilt = icond{i};
        % exclude burnin trials
        ifilt = ifilt(nburnin+1:end);
        nfilt = numel(ifilt);
        % get data
        sess = cat(1,sess,ceil(dat.block(ifilt)/8));
        strc = cat(1,strc,dat.staircase(ifilt));
        prop = cat(1,prop,dat.prop(ifilt));
        resp = cat(1,resp,dat.resp(ifilt));
        conf = cat(1,conf,dat.conf(ifilt));
    end
    
    % fit confidence model
    cfg      = [];
    cfg.stim = strc;
    cfg.resp = resp;
    cfg.conf = conf;
    cfg.sigcon = 0;
    cfg.nsmp = 1e2; % number of particle samples
    cfg.nres = 1e2; % number of bootstrap resamples
    out = fit_model_titration(cfg);
    out_fit{idata} = out;
    fprintf('done.\n');
end
% save output
save(fullfile(savepath,'expe2_results_fit_confidence_no_noise.mat'),'out_fit');

%% SUPPLEMENTARY: FIT CONFIDENCE MODEL TO TITRATION DATA - PER SESSION
fprintf('<strong>Fitting confidence model to titration data - per session</strong>\n\n')

clearvars -except datapath savepath ndata participants
% load titration results
load(fullfile(datapath,'expe2_titration.mat'));
dat = expe2_titration;
clear('expe2_titration');

% set analysis parameters
nburnin = 50; % number of burnin trials at the beginning of each session

% get list of participants
participants = unique(dat.participant);
ndata = numel(participants);

% initialize output
out_fit = cell(ndata,2); % fitting output

for idata = 1:ndata
    fprintf('\tParticipant %02d/%02d...',idata,ndata)
    % get data
    sess = []; % session number (1/2)
    strc = []; % titration side (+1:light/-1:dark)
    prop = []; % stimulus proportion of light
    resp = []; % response (+1:light/-1:dark)
    conf = []; % confidence (+1:high/-1:low)
    % define conditions
    icond = { ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase > 0 & dat.block > 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block <= 8), ...
        find(dat.participant == participants(idata) & ...
        dat.staircase < 0 & dat.block > 8)};
    ncond = numel(icond);
    for i = 1:ncond
        % get condition
        ifilt = icond{i};
        % exclude burnin trials
        ifilt = ifilt(nburnin+1:end);
        nfilt = numel(ifilt);
        % get data
        sess = cat(1,sess,ceil(dat.block(ifilt)/8));
        strc = cat(1,strc,dat.staircase(ifilt));
        prop = cat(1,prop,dat.prop(ifilt));
        resp = cat(1,resp,dat.resp(ifilt));
        conf = cat(1,conf,dat.conf(ifilt));
    end
    
    for isess = 1:2
        % filter current session
        ifilt = sess == isess;
        
        % fit confidence model
        cfg      = [];
        cfg.stim = strc(ifilt);
        cfg.resp = resp(ifilt);
        cfg.conf = conf(ifilt);
        cfg.nsmp = 1e2; % number of particle samples
        cfg.nres = 1e2; % number of bootstrap resamples
        out = fit_model_titration(cfg);
        out_fit{idata,isess} = out;
        
    end
    fprintf('done.\n');
end
% save output
save(fullfile(savepath,'expe2_results_fit_confidence_per_session.mat'),'out_fit');

%% SUPPLEMENTARY: FIT STABILIZING MODELS PER SESSION
fprintf('<strong>Fitting stabilizing models per session</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

models = {'senbias','inflaps','infdisc','selbias','selepsi'};
nmodel = numel(models);
out_fit = cell(ndata,nmodel,2); % fit session 1 and session 2 separately
for idata = 1:ndata
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    for isess = 1:2
        if isess == 1
            % filter participant
            ifilt = dat.participant == participants(idata) & dat.block <= 8;
        elseif isess == 2
            ifilt = dat.participant == participants(idata) & dat.block > 8;
        end
        
        % get data
        t = dat.trial(ifilt);   % trial number in block
        s = dat.stim(ifilt);    % stimulus (+1/-1)
        w = dat.wght(ifilt);    % strength (1/2/3)
        r = dat.resp(ifilt);    % response (+1/-1)
        
        % configure model fitting
        cfg_fit          = [];
        cfg_fit.t        = t;               % trial number in block
        cfg_fit.s        = s;               % stimulus
        cfg_fit.w        = w;               % strength
        cfg_fit.r        = r;               % response
        cfg_fit.nsmp     = 1e3;             % number of samples
        cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
        cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
        cfg_fit.sigsel   = 0;               % no selection noise
        
        for imodel = 1:nmodel
            % fit specified model
            fprintf('\tmodel %s...',models{imodel});
            cfg_tmp = cfg_fit;
            cfg_tmp.modtype = models{imodel};
            out_fit{idata,imodel,isess} = fit_model_revrepwgh(cfg_tmp);
        end
    end
end
% save results
save(fullfile(savepath,'expe2_results_fit_stabilizing_models_per_session.mat'),'out_fit');

%% SUPPLEMENTARY: FIT FLAT UNCONDITIONAL INFERENCE MODEL
fprintf('<strong>Fitting flat unconditional inference model</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    w = dat.wght(ifilt);    % strength (1/2/3)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;               % trial number in block
    cfg_fit.s        = s;               % stimulus
    cfg_fit.w        = w;               % strength
    cfg_fit.r        = r;               % response
    cfg_fit.nsmp     = 1e3;             % number of samples
    cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
    cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
    cfg_fit.sigsel   = 0;               % no selection noise
    
    % fit flat inference
    cfg_fit.modtype  = 'flatinf';
    out_fit{idata} = fit_model_revrepwgh(cfg_fit);
end
% save results
save(fullfile(savepath,'expe2_results_fit_flat_inference_model.mat'),'out_fit');

%% SUPPLEMENTARY: FIT HIERARCHICAL WEIGHTED INFERENCE MODEL
fprintf('<strong>Fitting weighted inference model</strong>\n')

clearvars -except datapath savepath ndata participants
load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
dat = expe2_data;
clear('expe2_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    w = dat.wght(ifilt);    % strength (1/2/3)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;               % trial number in block
    cfg_fit.s        = s;               % stimulus
    cfg_fit.w        = w;               % strength
    cfg_fit.r        = r;               % response
    cfg_fit.nsmp     = 1e3;             % number of samples
    cfg_fit.nres     = 1e3;             % number of re-samples for bootstrapping
    cfg_fit.pfal     = [0.3 0.2 0.1];   % false-alarm rate
    cfg_fit.sigsel   = 0;               % no selection noise
    cfg_fit.modtype  = 'infwght';
    
    out_fit{idata} = fit_model_revrepwgh(cfg_fit);
end
% save results
save(fullfile(savepath,'expe2_results_fit_weighted_inference_model.mat'),'out_fit');

