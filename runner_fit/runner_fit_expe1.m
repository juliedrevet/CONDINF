%% RUNNER FIT EXPERIMENT 1
% needed on the path for fitting models to participants' behaviour:
%       *VBMC v1.0 (https://github.com/lacerbi/vbmc/)
%       *BADS v1.0 (https://github.com/lacerbi/bads/)

% clear workspace
clearvars
close all
clc

addpath '../matlab_functions'; % all fitting functions are stored in this folder

datapath = '../data_table/';
savepath = '../results/';
if ~exist(savepath,'dir')
    mkdir(savepath);
end
%% check participants 
nstd = 3.5;
pcor = [];
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');
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

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

out_fit = cell(ndata,1);
for idata = 1:ndata    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    cfg_fit.modtype  = 'nonstab';
    
    out_fit{idata} = fit_model_revrep(cfg_fit);
end
% save results
save(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit');

%% FIT STABILIZING PARAMETERS BEST-FITTING PARTICIPANTS' OVERALL SWITCH RATE
fprintf('<strong>Fitting stabilizing parameters best-fitting participants'' overall switch rate</strong>\n')

clearvars -except ndata participants savepath datapath

% start from unstabilized best-fitting model parameters and add the
% stabilizing parameter
load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit')
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
        out_fit{idata,imodel} = fit_overall_stab(cfg_fit);
        fprintf('done.\n')
    end
end
% save results
save(fullfile(savepath,'expe1_results_fit_stabilizing_parameters.mat'),'out_fit');

%% EX-ANTE RECOVERY STABILIZING MODELS
fprintf('<strong>Ex-ante recovery between stabilizing models</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

models = {'senbias','inflaps','infdisc','selbias','selepsi'};
nmodel = numel(models);
nsim = 10; % number of simulations to be recovered

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
        t = dat.trial(ifilt);   % trial number in block
        s = dat.stim(ifilt);    % stimulus (+1/-1)
        for imodel = 1:nmodel
            fprintf('processing participant %02d/%02d - simulation model %s - simulation %02d\n\n',idata,ndata,models{imodel},isim);
            modtype = models{imodel};
            
            % 1/ simulate model
            cfg_sim = [];
            cfg_sim.modtype = modtype;
            cfg_sim.t = t;
            cfg_sim.s = s;
            cfg_sim.pfal = 0.20;
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
            out = sim_model(cfg_sim);
            
            % configure model fitting
            cfg_fit          = [];
            cfg_fit.t        = cfg_sim.t;   % trial number in block
            cfg_fit.s        = cfg_sim.s;   % stimulus
            cfg_fit.r        = out.r;       % response
            cfg_fit.nsmp     = 1e3;         % number of samples
            cfg_fit.nres     = 1e3;         % number of re-samples for bootstrapping
            cfg_fit.pfal     = 0.2;         % false-alarm rate
            cfg_fit.verbose  = false;       % verbose VBMC output
            cfg_fit.sigsel   = 0;
            
            for imodrec = 1:nmodel
                fprintf('\trecovery with model %s\n\n',models{imodrec});
                cfg_tmp = cfg_fit;
                cfg_tmp.modtype = models{imodrec};
                out_fit = fit_model_revrep(cfg_tmp);
                elbo(idata,imodel,imodrec,isim) = out_fit.elbo;
            end
        end
    end
end
% save results
save(fullfile(savepath,'expe1_results_recovery_stabilizing_models_elbo.mat'),'elbo');

%% FIT STABILIZING MODELS
fprintf('<strong>Fitting stabilizing models</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

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
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    
    for imodel = 1:nmodel
        % fit specified model
        fprintf('\tmodel %s...',models{imodel});
        cfg_tmp = cfg_fit;
        cfg_tmp.modtype = models{imodel};
        out_fit{idata,imodel} = fit_model_revrep(cfg_tmp);
    end
end
% save results
save(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');

%% SUPPLEMENTARY: EX-ANTE RECOVERY BASE MODEL SELECTION VERSUS INFERENCE NOISE
fprintf('<strong>Ex-ante recovery base model - selection versus inference noise</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

% output structure
out_rec = cell(2,2);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    
    % get simulation parameters from their a-priori distribution
    % 1/ perceived hazard rate
    h0 = betarnd(2,18);
    % 2/ inference noise
    siginf0 = gamrnd(4,0.25);
    % 3/ selection noise
    sigsel0 = gamrnd(4,0.25);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    
    % 1/ simulate (non-stabilized) INFERENCE noise model
    cfg_sim = [];
    cfg_sim.modtype = 'nonstab';
    cfg_sim.nsmp    = 1;
    cfg_sim.t       = t;        % trial number in block
    cfg_sim.s       = s;        % stimulus
    cfg_sim.pfal    = 0.2;      % false-alarm rate
    cfg_sim.h       = h0;       % perceived hazard rate
    cfg_sim.siginf  = siginf0;  % inference noise
    cfg_sim.sigsel  = 0;        % no selection noise
    cfg_sim.seed    = idata*100;
    out = sim_model(cfg_sim);
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.modtype  = 'nonstab';
    cfg_fit.t        = cfg_sim.t;   % trial number in block
    cfg_fit.s        = cfg_sim.s;   % stimulus
    cfg_fit.r        = out.r;       % response
    cfg_fit.pfal     = 0.2;         % false-alarm rate
    cfg_fit.nsmp     = 1e3;         % number of samples
    cfg_fit.nres     = 1e3;         % number of re-samples for bootstrapping
    cfg_fit.verbose  = false;       % verbose VBMC output
    
    % 1.1/ fit (non-stabilized) INFERENCE noise model
    cfg_tmp = cfg_fit;
    cfg_tmp.sigsel = 0; % no selection noise
    out_rec{1,1} = fit_model_revrep(cfg_tmp);
    
    % 1.2/ fit (non-stabilized) SELECTION noise model
    cfg_tmp = cfg_fit;
    cfg_tmp.siginf = 0; % no inference noise
    out_rec{1,2} = fit_model_revrep(cfg_tmp);
    
    % 2/ simulate (non-stabilized) SELECTION noise model
    cfg_sim = [];
    cfg_sim.modtype = 'nonstab';
    cfg_sim.nsmp    = 1;
    cfg_sim.t       = t;        % trial number in block
    cfg_sim.s       = s;        % stimulus
    cfg_sim.pfal    = 0.2;      % false-alarm rate
    cfg_sim.h       = h0;       % perceived hazard rate
    cfg_sim.siginf  = 0;        % no inference noise
    cfg_sim.sigsel  = sigsel0;  % selection noise
    cfg_sim.seed    = idata*100;
    out = sim_model(cfg_sim);
    
    % configure model fitting
    cfg_fit         = [];
    cfg_fit.modtype = 'nonstab';
    cfg_fit.t       = cfg_sim.t;    % trial number in block
    cfg_fit.s       = cfg_sim.s;    % stimulus
    cfg_fit.r       = out.r;        % response
    cfg_fit.pfal    = 0.2;          % false-alarm rate
    cfg_fit.nsmp    = 1e3;          % number of samples
    cfg_fit.nres    = 1e3;          % number of re-samples for bootstrapping
    cfg_fit.verbose = false;        % verbose VBMC output
    
    % 2.1/ fit (non-stabilized) INFERENCE noise model
    cfg_tmp = cfg_fit;
    cfg_tmp.sigsel = 0; % no selection noise
    out_rec{2,1} = fit_model_revrep(cfg_tmp);
    
    % 2.2/ fit (non-stabilized) SELECTION noise model
    cfg_tmp = cfg_fit;
    cfg_tmp.siginf = 0; % no inference noise
    out_rec{2,2} = fit_model_revrep(cfg_tmp);
end
% save results
save(fullfile(savepath,'expe1_results_recovery_infvsel_base_model.mat'),'out_rec');

%% SUPPLEMENTARY: FIT SELECTION NOISE MODELS
fprintf('<strong>Fitting noisy selection models</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

out_fit = cell(ndata,2);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.siginf   = 0;   % no inference noise
    
    % fit model without stabilization
    cfg_fit.modtype  = 'nonstab';
    out_fit{idata,1} = fit_model_revrep(cfg_fit);
    % fit conditional inference model
    cfg_fit.modtype  = 'infdisc';
    out_fit{idata,2} = fit_model_revrep(cfg_fit);   
end
% save results
save(fullfile(savepath,'expe1_results_fit_selection_noise_models.mat'),'out_fit');

%% SUPPLEMENTARY: EX-ANTE RECOVERY CONDITIONAL INFERENCE MODEL SELECTION VERSUS INFERENCE NOISE
fprintf('<strong>Ex-ante recovery conditional inference model - selection versus inference noise</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

% output structure
out_rec = cell(ndata,2,2);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d\n',idata,ndata);
    
    % get simulation parameters from their a-priori distribution
    % 1/ perceived hazard rate
    h0 = betarnd(2,18);
    % 2/ inference noise
    siginf0 = gamrnd(4,0.25);
    % 3/ selection noise
    sigsel0 = gamrnd(4,0.25);
    % 4/ reliability threshold delta
    delta0 = gamrnd(1,0.75);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    
    % 1/ simulate conditional inference model with INFERENCE noise
    cfg_sim = [];
    cfg_sim.modtype = 'infdisc';
    cfg_sim.nsmp    = 1;
    cfg_sim.t       = t;        % trial number in block
    cfg_sim.s       = s;        % stimulus
    cfg_sim.pfal    = 0.2;      % false-alarm rate
    cfg_sim.h       = h0;       % perceived hazard rate
    cfg_sim.siginf  = siginf0;  % inference noise
    cfg_sim.sigsel  = 0;        % selection noise
    cfg_sim.delta   = delta0;   % reliability threshold
    cfg_sim.seed    = idata*100;
    out = sim_model(cfg_sim);
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.modtype  = 'infdisc';
    cfg_fit.t        = cfg_sim.t;   % trial number in block
    cfg_fit.s        = cfg_sim.s;   % stimulus
    cfg_fit.r        = out.r;       % response
    cfg_fit.pfal     = 0.2;         % false-alarm rate
    cfg_fit.nsmp     = 1e3;         % number of samples
    cfg_fit.nres     = 1e3;         % number of re-samples for bootstrapping
    cfg_fit.verbose  = false;       % verbose VBMC output
    
    % 1.1/ fit conditional inference model with INFERENCE noise
    cfg_tmp = cfg_fit;
    cfg_tmp.sigsel = 0; % no selection noise
    out_rec{idata,1,1} = fit_model_revrep(cfg_tmp);
    
    % 1.2/ fit conditional inference model with SELECTION noise
    cfg_tmp = cfg_fit;
    cfg_tmp.siginf = 0; % no inference noise
    out_rec{idata,1,2} = fit_model_revrep(cfg_tmp);
    
    % 2/ simulate conditional inference model with SELECTION noise
    cfg_sim = [];
    cfg_sim.modtype = 'infdisc';
    cfg_sim.nsmp    = 1;
    cfg_sim.t       = t;        % trial number in block
    cfg_sim.s       = s;        % stimulus
    cfg_sim.pfal    = 0.2;      % false-alarm rate
    cfg_sim.h       = h0;       % preceived hazard rate
    cfg_sim.siginf  = 0;        % inference noise
    cfg_sim.sigsel  = sigsel0;  % selection noise
    cfg_sim.delta   = delta0;   % reliability threshold
    cfg_sim.seed    = idata*100;
    out = sim_model(cfg_sim);
    
    % configure model fitting
    cfg_fit         = [];
    cfg_fit.modtype = 'infdisc';
    cfg_fit.t       = cfg_sim.t;    % trial number in block
    cfg_fit.s       = cfg_sim.s;    % stimulus
    cfg_fit.r       = out.r;        % response
    cfg_fit.pfal    = 0.2;          % false-alarm rate
    cfg_fit.nsmp    = 1e3;          % number of samples
    cfg_fit.nres    = 1e3;          % number of re-samples for bootstrapping
    cfg_fit.verbose = false;        % verbose VBMC output
    
    % 2.1/ fit conditional inference model with INFERENCE noise
    cfg_tmp = cfg_fit;
    cfg_tmp.sigsel = 0; % no selection noise
    out_rec{idata,2,1} = fit_model_revrep(cfg_tmp);
    
    % 2.2/ fit conditional inference model with SELECTION noise
    cfg_tmp = cfg_fit;
    cfg_tmp.siginf = 0; % no inference noise
    out_rec{idata,2,2} = fit_model_revrep(cfg_tmp);
end
% save results
save(fullfile(savepath,'expe1_results_recovery_infvsel_conditional_inference_model.mat'),'out_rec');

%% SUPPLEMENTARY: FIT LEAKY INFERENCE MODELS
fprintf('<strong>Fitting leaky inference models</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

% leaky inference w/o stabilization
fprintf('no stabilization\n')
out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    cfg_fit.modtype  = 'nonstab';
    
    out_fit{idata} = fit_model_revrep_leak(cfg_fit);
end
% save results
save(fullfile(savepath,'expe1_results_fit_leak_noisy_inference_model.mat'),'out_fit');

% leaky inference with stabilization
fprintf('leaky inference stabilizing models\n')
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
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    
    for imodel = 1:nmodel
        % fit specified model
        fprintf('\tmodel %s...',models{imodel});
        cfg_tmp = cfg_fit;
        cfg_tmp.modtype = models{imodel};
        out_fit{idata,imodel} = fit_model_revrep_leak(cfg_tmp);
    end
end
% save results
save(fullfile(savepath,'expe1_results_fit_leak_stabilizing_models.mat'),'out_fit');

%% SUPPLEMENTARY: FIT STABILIZING MODELS PER SESSION
fprintf('<strong>Fitting stabilizing models per session</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

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
        r = dat.resp(ifilt);    % response (+1/-1)
        
        % configure model fitting
        cfg_fit          = [];
        cfg_fit.t        = t;   % trial number in block
        cfg_fit.s        = s;   % stimulus
        cfg_fit.r        = r;   % response
        cfg_fit.nsmp     = 1e3; % number of samples
        cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
        cfg_fit.pfal     = 0.2; % false-alarm rate
        cfg_fit.sigsel   = 0;   % no selection noise
        
        for imodel = 1:nmodel
            % fit specified model
            fprintf('\tmodel %s...',models{imodel});
            cfg_tmp = cfg_fit;
            cfg_tmp.modtype = models{imodel};
            out_fit{idata,imodel,isess} = fit_model_revrep(cfg_tmp);
        end
    end
end
% save results
save(fullfile(savepath,'expe1_results_fit_stabilizing_models_per_session.mat'),'out_fit');

%% SUPPLEMENTARY: FIT FLAT UNCONDITIONAL INFERENCE MODEL
fprintf('<strong>Fitting flat unconditional inference model</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    
    % fit flat inference
    cfg_fit.modtype  = 'flatinf';
    out_fit{idata} = fit_model_revrep(cfg_fit);
end
% save results
save(fullfile(savepath,'expe1_results_fit_flat_inference_model.mat'),'out_fit');

%% SUPPLEMENTARY: FIT HIERARCHICAL WEIGHTED INFERENCE MODEL
fprintf('<strong>Fitting weighted inference model</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

out_fit = cell(ndata,1);
for idata = 1:ndata
    
    fprintf('processing participant %02d/%02d...',idata,ndata);
    
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    t = dat.trial(ifilt);   % trial number in block
    s = dat.stim(ifilt);    % stimulus (+1/-1)
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit          = [];
    cfg_fit.t        = t;   % trial number in block
    cfg_fit.s        = s;   % stimulus
    cfg_fit.r        = r;   % response
    cfg_fit.nsmp     = 1e3; % number of samples
    cfg_fit.nres     = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal     = 0.2; % false-alarm rate
    cfg_fit.sigsel   = 0;   % no selection noise
    cfg_fit.modtype  = 'infwght';
    
    out_fit{idata} = fit_model_revrep(cfg_fit);
end
% save results
save(fullfile(savepath,'expe1_results_fit_weighted_inference_model.mat'),'out_fit');

%% SUPPLEMENTARY: FIT STABILIZING MODELS TO EACH CHOICE
fprintf('<strong>Fitting stabilizing models to each choice</strong>\n')

clearvars -except ndata participants savepath datapath
load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
dat = expe1_data;
clear('expe1_data');

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
    r = dat.resp(ifilt);    % response (+1/-1)
    
    % configure model fitting
    cfg_fit         = [];
    cfg_fit.t       = t;   % trial number in block
    cfg_fit.s       = s;   % stimulus
    cfg_fit.r       = r;   % response
    cfg_fit.nsmp    = 1e3; % number of samples
    cfg_fit.nres    = 1e3; % number of re-samples for bootstrapping
    cfg_fit.pfal    = 0.2; % false-alarm rate
    cfg_fit.sigsel  = 0;   % no selection noise
    
    for imodel = 1:nmodel
        % fit specified model
        fprintf('\tmodel %s...',models{imodel});
        cfg_tmp = cfg_fit;
        cfg_tmp.modtype = models{imodel};
        if ~strcmp(cfg_tmp.modtype,'selbias')
            cfg_tmp.beta = 0;
        end
        out_fit{idata,imodel} = fit_model_choice(cfg_tmp);
        out_fit{idata,imodel} = rmfield(out_fit{idata,imodel},{'p','x'}); % not necessary here
        % add reversal metrics (participant)
        c = getc(cfg_tmp.t,cfg_tmp.s,cfg_tmp.r);
        out_fit{idata,imodel}.crev_sub = c.rev_avg;
        out_fit{idata,imodel}.crep_sub = c.rep_avg;
        % add reversal metrics (model predictions)
        cfg_sim = cfg_tmp;
        cfg_sim.h = out_fit{idata,imodel}.h;
        cfg_sim.siginf = out_fit{idata,imodel}.siginf;
        cfg_sim.sigsel = out_fit{idata,imodel}.sigsel;
        cfg_sim.epsi = out_fit{idata,imodel}.epsi;
        cfg_sim.beta = out_fit{idata,imodel}.beta;
        cfg_sim.(out_fit{idata,imodel}.xnam{end}) = out_fit{idata,imodel}.(out_fit{idata,imodel}.xnam{end}); % add stab param
        out_sim = sim_model(cfg_sim);
        c = getc(cfg_tmp.t,cfg_tmp.s,out_sim.r);
        out_fit{idata,imodel}.crev_avg = c.rev_avg;
        out_fit{idata,imodel}.crev_std = c.rev_std;
        out_fit{idata,imodel}.crep_avg = c.rep_avg;
        out_fit{idata,imodel}.crep_std = c.rep_std;
    end
end
% save results
save(fullfile(savepath,'expe1_results_fit_stabilizing_models_each_choice.mat'),'out_fit');