 %% RUNNER ANALYSIS SUPPLEMENTARY FIGURES
% this script runs analysis + plots from supplementary
% needed on the path for running analysis:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)

addpath '../helper_plot';
addpath '../matlab_functions';

%% SUPPLEMENTARY: COMPARISON BETWEEN INFERENCE AND SELECTION NOISE (Supplementary Figure 1)
fprintf('<strong>*** SUPPLEMENTARY: COMPARISON BETWEEN INFERENCE AND SELECTION NOISE ***</strong>\n\n');

% BASE MODEL
clearvars;
savepath = '../results';
% ex-ante confusion matrix on simulated data
load(fullfile(savepath,'expe1_results_recovery_infvsel_base_model.mat'),'out_rec')
elbo_rec = cellfun(@(s)getfield(s,'elbo'),out_rec);
nmod = size(elbo_rec,2);

alphac = nan(nmod,nmod);
pexc   = nan(nmod,nmod);
for imodsim = 1:nmod
    alphac(imodsim,:) = spm_BMS(squeeze(elbo_rec(:,imodsim,:))); % get alpha from spm_bms for each simulation
end
% compute by hand exceedance probability from each originating simulation
% models
x_col = nan(nmod,nmod,1e6);
for imodrec = 1:nmod
    x_col(:,imodrec,:) = drchrnd(alphac(:,imodrec),1e6)';
    [~,J] = max(x_col(:,imodrec));
    for imodsim = 1:nmod
        pexc(imodsim,imodrec) = mean(J==imodsim);
    end
end

figname = 'expe1_recovery_infvsel_base_model';
pbar = 1; HGT = 4;
figure('Color','white','Name',figname);
cmap = gray; cmap = flipud(cmap);
imagesc(nanmean(pexc,3),[0 1]); hold on
colormap(cmap); colorbar('off');
ylabel('simulated model');xlabel('recovered model');
set(gca,'TickDir','ou','XTick',1:nmod,'YTick',1:nmod)
set(gca,'XTickLabel',{'inference','selection'})
set(gca,'YTickLabel',{'inference','selection'})
for k = 1:nmod
    text(k-.25,k,sprintf('>%3.3f',floor(pexc(k,k)*1000)/1000),'FontSize',8,'Color','w')
end
helper_ILL(gcf,pbar,HGT);
set(gca,'Layer','top','Box','on')

% BAYESIAN MODEL SELECTION on human data
load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit')
out_fit_infvsel = out_fit;
clearvars out_fit;
load(fullfile(savepath,'expe1_results_fit_selection_noise_models.mat'),'out_fit')
out_fit_infvsel = cat(2,out_fit_infvsel,out_fit(:,1));
elbo_infvsel = cellfun(@(s)getfield(s,'elbo'),out_fit_infvsel);

rgb = [0.50,0.50,0.50;... 
       0.75 0.75 0.75];

cfg = [];
cfg.figname = 'expe1_bms';
cfg.rgb = rgb;
cfg.elbo = elbo_infvsel;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
cfg.models = {'\sigma_{inf}' '\sigma_{sel}'};
helper_plot_bms(cfg);

% CONDITIONAL INFERENCE MODEL
clearvars;
savepath = '../results';
% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapse
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapse
    ]/255;
% ex-ante confusion matrix on simulated data
load(fullfile(savepath,'expe1_results_recovery_infvsel_conditional_inference_model.mat'),'out_rec')
elbo_rec = cellfun(@(s)getfield(s,'elbo'),out_rec);
nmod = size(elbo_rec,2);

alphac = nan(nmod,nmod);
pexc   = nan(nmod,nmod);
for imodsim = 1:nmod
    alphac(imodsim,:) = spm_BMS(squeeze(elbo_rec(:,imodsim,:))); % get alpha from spm_bms for each simulation
end
% compute by hand exceedance probability from each originating simulation
% models
x_col = nan(nmod,nmod,1e6);
for imodrec = 1:nmod
    x_col(:,imodrec,:) = drchrnd(alphac(:,imodrec),1e6)';
    [~,J] = max(x_col(:,imodrec));
    for imodsim = 1:nmod
        pexc(imodsim,imodrec) = mean(J==imodsim);
    end
end

figname = 'expe1_recovery_infvsel_conditional_inference_model';
pbar = 1; HGT = 4;
figure('Color','white','Name',figname);
cmap = gray; cmap = flipud(cmap);
imagesc(nanmean(pexc,3),[0 1]); hold on
colormap(cmap); colorbar('off');
ylabel('simulated model');xlabel('fitted model');
set(gca,'TickDir','ou','XTick',1:nmod,'YTick',1:nmod)
set(gca,'XTickLabel',{'inference','selection'})
set(gca,'YTickLabel',{'inference','selection'})
for k = 1:nmod
    text(k-.25,k,sprintf('>%3.3f',floor(pexc(k,k)*1000)/1000),'FontSize',8,'Color','w')
end
fig = helper_ILL(gcf,pbar,HGT);
set(gca,'Layer','top','Box','on')

% BAYESIAN MODEL SELECTION on human data
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit')
out_fit_infvsel = out_fit(:,3); % conditional inference model
clearvars out_fit;
load(fullfile(savepath,'expe1_results_fit_selection_noise_models.mat'),'out_fit')
out_fit_infvsel = cat(2,out_fit_infvsel,out_fit(:,2));
elbo_infvsel = cellfun(@(s)getfield(s,'elbo'),out_fit_infvsel);

rgb_infvsel = rgb(3,:);
rgb_infvsel(2,:) = (rgb_infvsel+1)/2;

cfg = [];
cfg.figname = 'expe1_bms';
cfg.rgb = rgb_infvsel;
cfg.elbo = elbo_infvsel;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
cfg.models = {'\sigma_{inf}' '\sigma_{sel}'};
helper_plot_bms(cfg);

%% SUPPLEMENTARY: COMPARISON BETWEEN LEAKY AND NORMATIVE EVIDENCE ACCUMULATION (Supplementary Figure 2)
fprintf('<strong>*** SUPPLEMENTARY: COMPARISON BETWEEN LEAKY AND NORMATIVE EVIDENCE ACCUMULATION ***</strong>\n\n');

clearvars;
savepath = '../results';

% Supplementary Figure 2a
fprintf('Comparing leaky and normative integration models (non-stabilized models)\n')
% define model color for plots
rgb = [0.5 0.5 0.5]; 

% extract model fitting output (without stabilization)
load(fullfile(savepath,'expe1_results_fit_leak_noisy_inference_model.mat'),'out_fit');
out_fit_leak = out_fit; clearvars out_fit;
load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit');
out_fit_norm = out_fit; clearvars out_fit;

out_fit = cat(2,out_fit_leak,out_fit_norm);
ndata = size(out_fit,1);
nmodel = size(out_fit,2);
modnam = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);

xavg = cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',0);
xstd = cellfun(@(x)getfield(x,'xstd'),out_fit,'UniformOutput',0);
xnam = out_fit{1}.xnam;
npar = numel(xnam);

elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_mod = cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0);
crep_mod = cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0);
% participants
crev_dat = cell2mat(cellfun(@(s)getfield(s,'crev_sub'),out_fit(:,1),'UniformOutput',0)')';
crep_dat = cell2mat(cellfun(@(s)getfield(s,'crep_sub'),out_fit(:,1),'UniformOutput',0)')';

rgb_mod = [(rgb+1)/2; rgb];

cfg = [];
cfg.figname = 'expe1_leak_bms_nonstab';
cfg.rgb = rgb_mod;
cfg.elbo = elbo;
cfg.ylabel = 'model probability';
cfg.xlabel = 'accumulation';
cfg.models = {'leaky' 'normative'};
helper_plot_bms(cfg);

for imodel = 1:nmodel
    % plot response switch curve  
    cfg = [];
    cfg.figname = sprintf('expe1_leak_SWI_model_%s_%d',modnam{imodel},imodel);
    cfg.rgb = rgb_mod(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-mean(cell2mat(crep_mod(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crep_mod(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = 1-mean(crep_dat)';
    cfg.Yerr = std(crep_dat,[],1)'/sqrt(ndata);
    helper_plot_rep(cfg)
    % plot response reversal curve 
    cfg = [];
    cfg.figname = sprintf('expe1_leak_REV_model_%s_%d',modnam{imodel},imodel);
    cfg.rgb = rgb_mod(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = mean(cell2mat(crev_mod(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crev_mod(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = mean(crev_dat)';
    cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
    helper_plot_rev(cfg)
end

for ipar = 1:npar
    X = cellfun(@(x)x(ipar),xavg(:,1)); % leaky
    Y = cellfun(@(x)x(ipar),xavg(:,2)); % normative
    Xstd = cellfun(@(x)x(ipar),xstd(:,1)); % leaky
    Ystd = cellfun(@(x)x(ipar),xstd(:,2)); % normative
    fprintf('\t%s',xnam{ipar});
    fprintf('\t%4.3f +/- %4.3f (leaky) ; \t%4.3f +/- %4.3f (normative)\n',mean(X),std(X,[],1)/sqrt(ndata),mean(Y),std(Y,[],1)/sqrt(ndata));
    
    cfg = [];
    cfg.figname = sprintf('expe1_leak_norm_corr_nonstab_%s',xnam{ipar});
    cfg.xavg = X; cfg.xstd = Xstd;
    cfg.yavg = Y; cfg.ystd = Ystd;
    cfg.rgb = rgb;
    
    if strcmp(xnam{ipar},'h')
        cfg.xlim = [0 .4]; cfg.xtick = 0:.1:.4; cfg.xticklabel = {'0' '0.1' '0.2' '0.3' '0.4'};
    elseif strcmp(xnam{ipar},'siginf')
        cfg.xlim = [0 4.5]; cfg.xtick = 0:1:4.5; cfg.xticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(xnam{ipar},'delta')
        cfg.xlim = [0 3.5]; cfg.xtick = 0:1.0:3.5; cfg.xticklabel = {'0' '1.0' '2.0' '3.0'};
    end
    cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;

    cfg.xlabel = sprintf('%s - leaky',xnam{ipar});
    cfg.ylabel = sprintf('%s - normative',xnam{ipar});
    helper_plot_corr(cfg);
end

% Supplementary Figure 2b
fprintf('Comparing leaky and normative integration models (stabilized through conditional inference):\n')
% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapses
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapses
    ]/255;
imodel = 3; % infdisc stabilization model
rgb = rgb(imodel,:); % define selected model color

load(fullfile(savepath,'expe1_results_fit_leak_stabilizing_models.mat'),'out_fit');
out_fit_leak = out_fit(:,imodel); clearvars out_fit;
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');
out_fit_norm = out_fit(:,imodel); clearvars out_fit;
out_fit = cat(2,out_fit_leak,out_fit_norm);
ndata = size(out_fit,1);
nmodel = size(out_fit,2);
modnam = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);

elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_mod = cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0);
crep_mod = cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0);
% parameter values
xavg = cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',0);
xstd = cellfun(@(x)getfield(x,'xstd'),out_fit,'UniformOutput',0);
xnam = out_fit{1}.xnam;
npar = numel(xnam);
% participants
crev_dat = cell2mat(cellfun(@(s)getfield(s,'crev_sub'),out_fit(:,1),'UniformOutput',0)')';
crep_dat = cell2mat(cellfun(@(s)getfield(s,'crep_sub'),out_fit(:,1),'UniformOutput',0)')';

rgb_mod = [(rgb+1)/2; rgb];

cfg = [];
cfg.figname = 'expe1_leak_bms_infdisc';
cfg.rgb = rgb_mod;
cfg.elbo = elbo;
cfg.ylabel = 'model probability';
cfg.models = {'leaky' 'normative'};
cfg.xlabel = 'accumulation';
helper_plot_bms(cfg);

for imodel = 1:nmodel % loop over accumulation type (leaky or normative)
    % plot stabilization models versus non stabilized model    
    cfg = [];
    cfg.figname = sprintf('expe1_leak_SWI_model_%s_%d',modnam{imodel},imodel);
    cfg.rgb = rgb_mod(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-mean(cell2mat(crep_mod(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crep_mod(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = 1-mean(crep_dat)';
    cfg.Yerr = std(crep_dat,[],1)'/sqrt(ndata);
    helper_plot_rep(cfg)
    
    cfg = [];
    cfg.figname = sprintf('expe1_leak_REV_model_%s_%d',modnam{imodel},imodel);
    cfg.rgb = rgb_mod(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = mean(cell2mat(crev_mod(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crev_mod(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = mean(crev_dat)';
    cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
    helper_plot_rev(cfg)
end

for ipar = 1:npar
    X = cellfun(@(x)x(ipar),xavg(:,1)); % leaky
    Y = cellfun(@(x)x(ipar),xavg(:,2)); % normative
    Xstd = cellfun(@(x)x(ipar),xstd(:,1)); % leaky
    Ystd = cellfun(@(x)x(ipar),xstd(:,2)); % normative
    fprintf('\t%s',xnam{ipar});
    fprintf('\t%4.3f +/- %4.3f (leaky) ; \t%4.3f +/- %4.3f (normative)\n',mean(X),std(X,[],1)/sqrt(ndata),mean(Y),std(Y,[],1)/sqrt(ndata));
    
    cfg = [];
    cfg.figname = sprintf('expe1_leak_norm_corr_infdisc_%s',xnam{ipar});
    cfg.xavg = X; cfg.xstd = Xstd;
    cfg.yavg = Y; cfg.ystd = Ystd;
    cfg.rgb = rgb;
    
    if strcmp(xnam{ipar},'h')
        cfg.xlim = [0 .4]; cfg.xtick = 0:.1:.4; cfg.xticklabel = {'0' '0.1' '0.2' '0.3' '0.4'};
    elseif strcmp(xnam{ipar},'siginf')
        cfg.xlim = [0 4.5]; cfg.xtick = 0:1:4.5; cfg.xticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(xnam{ipar},'delta')
        cfg.xlim = [0 3.5]; cfg.xtick = 0:1.0:3.5; cfg.xticklabel = {'0' '1.0' '2.0' '3.0'};
    end
    cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;

    cfg.xlabel = sprintf('%s - leaky',xnam{ipar});
    cfg.ylabel = sprintf('%s - normative',xnam{ipar});
    helper_plot_corr(cfg);
end

%% SUPPLEMENTARY: PREDICTED EFFECTS OF RESPONSE STABILIZATION STRATEGIES EXPERIMENT 2 (Supplementary Figure 3)
fprintf('<strong>*** SUPPLEMENTARY: PREDICTED EFFECTS OF RESPONSE STABILIZATION STRATEGIES EXPERIMENT 2 ***</strong>\n\n');

clearvars;
savepath = '../results';

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapses
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapses
    ]/255;

if exist(fullfile(savepath,'expe2_predicted_effects_stabilization.mat'),'file')
    load(fullfile(savepath,'expe2_predicted_effects_stabilization.mat'))
else
    load(fullfile(savepath,'expe2_results_fit_stabilizing_parameters.mat'))
    ndata = size(out_fit,1);
    modnam = cellfun(@(s)getfield(s,'modtype'),out_fit(1,:),'UniformOutput',0);
    nmodel = numel(modnam);
    parnam = {'lambda','plaps','delta','beta','epsi'};
    models = {'Perceptual bias','Inference lapses','Conditional inference','Repetition bias','Response lapses'};
    
    stab_params = nan(ndata,nmodel);
    for imodel = 1:nmodel
        stab_params(:,imodel) = cellfun(@(s)getfield(s,parnam{imodel}),out_fit(:,imodel));
    end
    clearvars out_fit;
    
    fprintf('<strong>Stabilization parameter values best-fitting participants'' overall switch rate:</strong>\n');
    for imodel = 1:nmodel
        fprintf('\t%s \t%s\t= %5.3f +/- %5.3f\n' ,models{imodel},parnam{imodel},mean(stab_params(:,imodel)),std(stab_params(:,imodel))/sqrt(ndata));
    end
    % Compute simulated effects of the 5 candidate strategies on reversal
    % behavior
    fprintf('<strong>Computing ex-ante stabilization model simulations</strong>\n')
    % extract model fitting output (without stabilization) to be stabilized
    load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit');
    % response switch curves and response reversal curves
    crev_w_stab_ex_ante = cell(ndata,nmodel);
    crep_w_stab_ex_ante = cell(ndata,nmodel);
    for imodel = 1:nmodel
        fprintf('\t%s \t..',models{imodel});
        for idata = 1:ndata
            fprintf('.')
            cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
            cfg_sim = [];
            cfg_sim.t = cfg.t;
            cfg_sim.s = cfg.s;
            cfg_sim.c = cfg.c;
            cfg_sim.w = cfg.w;
            cfg_sim.h = out_fit{idata}.h;
            cfg_sim.siginf = out_fit{idata}.siginf;
            cfg_sim.modtype = modnam{imodel};
            cfg_sim.sigsel = 0;
            cfg_sim.pfal = [0.30 0.20 0.10];
            cfg_sim.nsmp = 1e3;
            % add stabilizing parameter
            cfg_sim.(parnam{imodel}) = stab_params(idata,imodel);
            out_sim = sim_model_expe2(cfg_sim);
            crev_w_stab_ex_ante{idata,imodel} = out_sim.c.rev_w_avg;
            crep_w_stab_ex_ante{idata,imodel} = out_sim.c.rep_w_avg;
        end
        fprintf('done.\n')
    end
    fprintf('\n');
    % Compute simulated effects of the 5 candidate strategies on response 
    % switches just before and after a reversal
    npar_sim = 20; % number of sample of stabilizing parameter to simulate
    stab_param_max = [0.5 .7 2 1.5 .7];   % max stab parameter values for simulation
    stab_param_vec = nan(nmodel,npar_sim); % parameter values for simulation
    last_before_w_avg = nan(nmodel,npar_sim,3);
    first_after_w_avg = nan(nmodel,npar_sim,3);
    last_before_w_err = nan(nmodel,npar_sim,3);
    first_after_w_err = nan(nmodel,npar_sim,3);
    last_before_w_ex_ante = nan(ndata,nmodel,3);
    first_after_w_ex_ante = nan(ndata,nmodel,3);
    fprintf('<strong>Computing fraction switch on last trial before and first trial after each reversal (long)</strong>\n')
    for imodel = 1:nmodel
        stab_param_vec(imodel,:) = linspace(0,stab_param_max(imodel),npar_sim);
        fprintf('\t%s \t..',models{imodel});
        for isim = 1:npar_sim
            fprintf('.')
            last_before_tmp = nan(ndata,3);
            first_after_tmp = nan(ndata,3);
            for idata = 1:ndata
                cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
                cfg_sim = [];
                cfg_sim.t = cfg.t;
                cfg_sim.s = cfg.s;
                cfg_sim.c = cfg.c;
                cfg_sim.w = cfg.w;
                cfg_sim.h = out_fit{idata}.h;
                cfg_sim.siginf = out_fit{idata}.siginf;
                cfg_sim.modtype = modnam{imodel};
                cfg_sim.sigsel = 0;
                cfg_sim.pfal = [0.30 0.20 0.10];
                cfg_sim.nsmp = 1e3;
                % add stabilizing parameter
                cfg_sim.(parnam{imodel}) = stab_param_vec(imodel,isim);
                out_sim = sim_model_expe2(cfg_sim);
                % extract fraction switch
                last_before_tmp(idata,:) = 1-out_sim.c.rep_w_avg(4,:);
                first_after_tmp(idata,:) = 1-out_sim.c.rep_w_avg(5,:);
            end
            last_before_w_avg(imodel,isim,:) = mean(last_before_tmp,1);
            first_after_w_avg(imodel,isim,:) = mean(first_after_tmp,1);
            last_before_w_err(imodel,isim,:) = std(last_before_tmp,[],1)/sqrt(ndata);
            first_after_w_err(imodel,isim,:) = std(first_after_tmp,[],1)/sqrt(ndata);
        end
        fprintf('done\n')
        % predictions for mean stabilization parameters matching overall
        % accuracy
        for idata = 1:ndata
            cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
            cfg_sim = [];
            cfg_sim.t = cfg.t;
            cfg_sim.s = cfg.s;
            cfg_sim.c = cfg.c;
            cfg_sim.w = cfg.w;
            cfg_sim.h = out_fit{idata}.h;
            cfg_sim.siginf = out_fit{idata}.siginf;
            cfg_sim.modtype = modnam{imodel};
            cfg_sim.sigsel = 0;
            cfg_sim.pfal = [0.30 0.20 0.10];
            cfg_sim.nsmp = 1e3;
            % add stabilizing mparameter
            cfg_sim.(parnam{imodel}) = mean(stab_params(:,imodel));
            out_sim = sim_model_expe2(cfg_sim);
            % extract fraction switch
            last_before_w_ex_ante(idata,imodel,:) = 1-out_sim.c.rep_w_avg(4,:);
            first_after_w_ex_ante(idata,imodel,:) = 1-out_sim.c.rep_w_avg(5,:);
        end
    end
    fprintf('\n');
    save(fullfile(savepath,'expe2_predicted_effects_stabilization.mat'),...
        'ndata','modnam','parnam','nmodel','models','stab_params',...
        'crev_w_stab_ex_ante','crep_w_stab_ex_ante','last_before_w_ex_ante','first_after_w_ex_ante',...
        'stab_param_vec',...
        'last_before_w_avg','first_after_w_avg','last_before_w_err','first_after_w_err')
end

% extract model fitting output (without stabilization) to be stabilized
load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit');
crev_w_nonstab = cellfun(@(s)getfield(s,'crev_w_avg'),out_fit,'UniformOutput',0);
crep_w_nonstab = cellfun(@(s)getfield(s,'crep_w_avg'),out_fit,'UniformOutput',0);

% Plot simulated effects of the 5 candidate strategies on reversal behavior
fprintf('<strong>Stabilization parameter values best-fitting participants'' overall switch rate:</strong>\n');
for imodel = 1:nmodel
    rgb_mod0 = rgb(imodel,:);
    rgb_dashed0 = [.5 .5 .5];
    rgb_mod = nan(3,3);
    rgb_dashed = nan(3,3);
    for istr = 1:3
        rgb_mod(istr,:) = (rgb_mod0+(3-istr))/(4-istr);
        rgb_dashed(istr,:) = (rgb_dashed0+(3-istr))/(4-istr);
    end
    cfg = [];
    cfg.figname = sprintf('expe2_ex_ante_model_%s_SWI',modnam{imodel});
    cfg.rgb = rgb_mod;
    cfg.rgb_dashed = rgb_dashed;
    cfg.yavg = 1-[mean(cell2mat(cellfun(@(x)x(:,1),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2)];
    cfg.yerr = [std(cell2mat(cellfun(@(x)x(:,1),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crep_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2)]/sqrt(ndata);
    cfg.yavg_dashed = 1-[mean(cell2mat(cellfun(@(x)x(:,1),crep_w_nonstab','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crep_w_nonstab','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crep_w_nonstab','UniformOutput',0)),2)];
    helper_plot_rep(cfg)
    
    cfg = [];
    cfg.figname = sprintf('expe2_ex_ante_model_%s_REV',modnam{imodel});
    cfg.rgb = rgb_mod;
    cfg.rgb_dashed = rgb_dashed;
    cfg.yavg = [mean(cell2mat(cellfun(@(x)x(:,1),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),2)];
    cfg.yerr = [std(cell2mat(cellfun(@(x)x(:,1),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crev_w_stab_ex_ante(:,imodel)','UniformOutput',0)),[],2)]/sqrt(ndata);
    cfg.yavg_dashed = [mean(cell2mat(cellfun(@(x)x(:,1),crev_w_nonstab','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crev_w_nonstab','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crev_w_nonstab','UniformOutput',0)),2)];
    helper_plot_rev(cfg)

    fprintf('\t%s \t%s\t= %5.3f +/- %5.3f\n' ,models{imodel},parnam{imodel},mean(stab_params(:,imodel)),std(stab_params(:,imodel))/sqrt(ndata));
end
fprintf('\n')

% Plot simulated effects of the 5 candidate strategies on response switches just before and after a reversal
pbar = 1.1; HGT = 4;
for imodel = 1:nmodel
    figname = sprintf('expe2_ex_ante_model_%s_firstlast',modnam{imodel});
    figure('Color','white','Name',figname);
    xvec = stab_param_vec(imodel,:);
    xl = [min(xvec) max(xvec)];
    hold on
    xlim(xl);
    ylim([0,.6]);
    for istr = 1:3
        rgb_str = (rgb(imodel,:)+(3-istr))/(4-istr);
        Y1 = last_before_w_avg(imodel,:,istr);
        Y2 = first_after_w_avg(imodel,:,istr);
        Y1err = last_before_w_err(imodel,:,istr);
        Y2err = first_after_w_err(imodel,:,istr);
        patch([xvec,fliplr(xvec)],[Y1+Y1err,fliplr(Y1-Y1err)], ...
            0.5*(rgb_str+1),'EdgeColor','none');
        patch([xvec,fliplr(xvec)],[Y2+Y2err,fliplr(Y2-Y2err)], ...
            0.5*(rgb_str+1),'EdgeColor','none');
    end
    plot([mean(stab_params(:,imodel)),mean(stab_params(:,imodel))],ylim,'Color',rgb(imodel,:));
    for istr = 1:3
        rgb_str = (rgb(imodel,:)+(3-istr))/(4-istr);
        Y1 = last_before_w_avg(imodel,:,istr);
        Y2 = first_after_w_avg(imodel,:,istr);
        plot(xvec,Y1,'Color',rgb_str,'LineWidth',1,'LineStyle',':');
        plot(xvec,Y2,'Color',rgb_str,'LineWidth',1);
        plot(mean(stab_params(:,imodel)),mean(last_before_w_ex_ante(:,imodel,istr)),'ko','MarkerSize',4,'MarkerFaceColor',rgb_str,'MarkerEdgeColor','non');
        plot(mean(stab_params(:,imodel)),mean(first_after_w_ex_ante(:,imodel,istr)),'ko','MarkerSize',4,'MarkerFaceColor',rgb_str,'MarkerEdgeColor','non');
    end
    set(gca,'YTick',0:0.2:.6)%,'XTick',xt,'XTickLabel',XTL);
    xlabel(sprintf('%s',parnam{imodel})); 
    ylabel('fraction switch');
    fig = helper_ILL(gcf,pbar,HGT,4);
end

%% SUPPLEMENTARY: BELIEF STABILIZATION THROUGH CONDITIONAL INFERENCE EXPERIMENT 2 (Supplementary Figure 4)
fprintf('<strong>*** BELIEF STABILIZATION THROUGH CONDITIONAL INFERENCE EXPERIMENT 2 ***</strong>\n\n');
clearvars;
savepath = '../results';

% CONFUSION MATRIX for EX-ANTE MODEL RECOVERY 
load(fullfile(savepath,'expe2_results_recovery_stabilizing_models_elbo.mat'),'elbo')
nmod = size(elbo,2);

alphac = nan(nmod,nmod);    % fitted model probabilities per simulated model
pexc = nan(nmod,nmod);      % exceedance probability of simulated model per fitting model
mean_elbo = mean(elbo,4);   % simulation + recovery done 10 times, take mean elbo
nsmp = 1e6;                 % number of samples from Dirichlet distribution to compute exceedance probability

% for each simulated model, compute probability of candidate fitting models
for imodsim = 1:nmod 
    alphac(imodsim,:) = spm_BMS(squeeze(mean_elbo(:,imodsim,:))); 
end
% compute for each fitted (recovered) model the exceedance probability of each
% simulated model
x_col = nan(nmod,nmod,nsmp);
for imodrec = 1:nmod
    x_col(:,imodrec,:) = drchrnd(alphac(:,imodrec),1e6)';
    [~,J] = max(x_col(:,imodrec,:));
    for imodsim = 1:nmod
        pexc(imodsim,imodrec) = mean(J==imodsim);
    end
end
HGT = 4.5; pbar = 1;
figure('Color','white','Name','expe2_confusion_matrix_stabilizing_models');
cmap = gray; cmap = flipud(cmap);
imagesc(pexc,[0 1])
colormap(cmap); cb=colorbar; title(cb,'model frequency')
ylabel('simulated model');xlabel('fitted model');
for i = 1:nmod
    pexc_diag = nanmean(pexc(i,i,:),3);
    text(i-0.5,i,sprintf('>%5.3f',floor(1000*pexc_diag)/1000),'Color','w')
end
set(gca,'TickDir','ou','XTick',1:nmod,'YTick',1:nmod)
helper_ILL(gcf,pbar,HGT);
set(gca,'Layer','top','Box','on')

% STABILIZING MODELS - EXTRACT FITTING OUTPUT
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'));
models = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);
ndata = size(out_fit,1);
nmodel = size(out_fit,2);

elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_w_stab = cellfun(@(s)getfield(s,'crev_w_avg'),out_fit,'UniformOutput',0);
crep_w_stab = cellfun(@(s)getfield(s,'crep_w_avg'),out_fit,'UniformOutput',0);

% stabilization parameters value
xavg = cellfun(@(x)x(3),cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',0));
xnam = cellfun(@(x)x(3),cellfun(@(x)getfield(x,'xnam'),out_fit(1,:),'UniformOutput',0));

% extract model fitting output (without stabilization)
load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'));
crev_w_nonstab = cellfun(@(s)getfield(s,'crev_w_avg'),out_fit,'UniformOutput',0);
crep_w_nonstab = cellfun(@(s)getfield(s,'crep_w_avg'),out_fit,'UniformOutput',0);

% participants
crev_w_dat = cellfun(@(s)getfield(s,'crev_w_sub'),out_fit,'UniformOutput',0);
crep_w_dat = cellfun(@(s)getfield(s,'crep_w_sub'),out_fit,'UniformOutput',0);

% BAYESIAN MODEL SELECTION 
rgb = [ ...
    238,197,141; ... % senbias
    048,102,158; ... % inflaps
    134,198,226; ... % infdisc
    171,131,171; ... % selbias
    099,074,135  ... % selepsi
    ]/255;

cfg = [];
cfg.figname = 'expe2_bms';
cfg.rgb = rgb;
cfg.elbo = elbo;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
helper_plot_bms(cfg);


% REVERSAL BEHAVIOR (MODEL SIMULATIONS)
fprintf('<strong>Stabilization parameters</strong>\n')
for imodel = 1:nmodel
    rgb_mod0 = rgb(imodel,:);
    rgb_mod = nan(3,3);
    for istr = 1:3
        rgb_mod(istr,:) = (rgb_mod0+(3-istr))/(4-istr);
    end
    % plot stabilization models for each stimulus strength
    cfg = [];
    cfg.figname = sprintf('expe2_SWI_model_%s',models{imodel});
    cfg.rgb = rgb_mod;
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-[mean(cell2mat(cellfun(@(x)x(:,1),crep_w_stab(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crep_w_stab(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crep_w_stab(:,imodel)','UniformOutput',0)),2)];
    cfg.yerr = [std(cell2mat(cellfun(@(x)x(:,1),crep_w_stab(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crep_w_stab(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crep_w_stab(:,imodel)','UniformOutput',0)),[],2)]/sqrt(ndata);
    cfg.Yavg = 1-[mean(cell2mat(cellfun(@(x)x(:,1),crep_w_dat','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crep_w_dat','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crep_w_dat','UniformOutput',0)),2)];
    cfg.Yerr = [std(cell2mat(cellfun(@(x)x(:,1),crep_w_dat','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crep_w_dat','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crep_w_dat','UniformOutput',0)),[],2)]/sqrt(ndata);
    cfg.ylim = [0 .6];
    helper_plot_rep(cfg)
    
    cfg = [];
    cfg.figname = sprintf('expe2_REV_model_%s',models{imodel});
    cfg.rgb = rgb_mod;
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = [mean(cell2mat(cellfun(@(x)x(:,1),crev_w_stab(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crev_w_stab(:,imodel)','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crev_w_stab(:,imodel)','UniformOutput',0)),2)];
    cfg.yerr = [std(cell2mat(cellfun(@(x)x(:,1),crev_w_stab(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crev_w_stab(:,imodel)','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crev_w_stab(:,imodel)','UniformOutput',0)),[],2)]/sqrt(ndata);
    cfg.Yavg = [mean(cell2mat(cellfun(@(x)x(:,1),crev_w_dat','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,2),crev_w_dat','UniformOutput',0)),2),...
                mean(cell2mat(cellfun(@(x)x(:,3),crev_w_dat','UniformOutput',0)),2)];
    cfg.Yerr = [std(cell2mat(cellfun(@(x)x(:,1),crev_w_dat','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,2),crev_w_dat','UniformOutput',0)),[],2),...
                std(cell2mat(cellfun(@(x)x(:,3),crev_w_dat','UniformOutput',0)),[],2)]/sqrt(ndata);
    helper_plot_rev(cfg)

    fprintf('\t%s',xnam{imodel});
    fprintf('\t%4.3f +/- %4.3f\n',mean(xavg(:,imodel)),std(xavg(:,imodel),[],1)/sqrt(ndata));
end
fprintf('\n')

%% SUPPLEMENTARY: DECISION CONFIDENCE MODEL (Supplementary Figure 5)
fprintf('<strong>*** SUPPLEMENTARY: DECISION CONFIDENCE MODEL ***</strong>\n\n');

clearvars;
savepath = '../results';

load(fullfile(savepath,'expe2_results_fit_confidence_per_session.mat'));
out_fit_session = out_fit;
load(fullfile(savepath,'expe2_results_fit_confidence_no_noise.mat'));
out_fit_no_noise = out_fit;
load(fullfile(savepath,'expe2_results_fit_confidence.mat'),'out_fit');

ndata = size(out_fit,1);

nsmp = 1e4;
pconf_dat       = nan(ndata,1);   % mean confidence for participant data
pdiff_dat_avg   = nan(ndata,1);   % mean accuracy difference between confident and unconfident trials for participant data
pdiff_dat_std   = nan(ndata,1);   % estimated accuracy difference s.d. between confident and unconfident trials for participant data
pconf_mod       = nan(ndata,2);   % mean confidence for model predictions
pdiff_mod_avg   = nan(ndata,2);   % mean accuracy difference between confident and unconfident trials for models predictions
pdiff_mod_std   = nan(ndata,2);   % estimated accuracy difference s.d. between confident and unconfident trials for models predictions
pconf_sess      = nan(ndata,2);   % mean confidence for participant data per session
pdiff_sess_avg  = nan(ndata,2);   % mean accuracy difference between confident and unconfident trials for participant data per session
pdiff_sess_std  = nan(ndata,2);   % estimated accuracy difference s.d. between confident and unconfident trials for participant data per session

for idata = 1:ndata
    fprintf('Participant %02d/%02d...',idata,ndata)
    resp = out_fit{idata}.cfg.resp;
    stim = out_fit{idata}.cfg.stim;
    conf = out_fit{idata}.cfg.conf;
    % get best-fitting parameters
    sigsen = out_fit{idata}.sigsen;
    sigcon = out_fit{idata}.sigcon;
    thrcon = out_fit{idata}.thrcon;
    
    % compute response/confidence statistics for participant data
    pconf_dat(idata,1) = mean(conf > 0);
    paccu_conf_tmp = [mean(resp(conf < 0) == stim(conf < 0));...  % unconfident
                      mean(resp(conf > 0) == stim(conf > 0))];    % confident
    pdiff_dat_avg(idata) = diff(paccu_conf_tmp);
    pdiff_dat_std(idata) = sqrt( ...
        paccu_conf_tmp(1)*(1-paccu_conf_tmp(1))/nnz(conf < 0)+ ... % estimated s.d. depending on number of trials
        paccu_conf_tmp(2)*(1-paccu_conf_tmp(2))/nnz(conf > 0));
    
    % compute response/confidence statistics for model predictions:
    % complete model
    pconf_mod(idata,1) = mean(out_fit{idata}.pconf);
    % simulate model response/confidence
    s = stim+sigsen*randn(length(stim),nsmp);                   % stimulus sensory response according to sensory noise
    r = sign(s);                                                % model response (+1/-1)
    c = sign(abs(s)+sigcon*randn(length(stim),nsmp)-thrcon);    % model confidence (+1/-1)
    
    stim_mat = repmat(stim,1,nsmp);
    paccu_conf_tmp = [mean(r(c<0) == stim_mat(c<0));...
                      mean(r(c>0) == stim_mat(c>0))];  
    pdiff_mod_avg(idata,1) = diff(paccu_conf_tmp);
    paccu_conf_tmp = ([sum((r==stim_mat)&c>0)./sum(c>0);...
                          sum((r==stim_mat)&c<0)./sum(c<0)]);
    pdiff_mod_std(idata,1) = sqrt(sum(var(paccu_conf_tmp,0,2)));
    
    % model without confidence noise
    pconf_mod(idata,2) = mean(out_fit_no_noise{idata}.pconf);
    % get best-fitting parameters
    sigsen = out_fit_no_noise{idata}.sigsen;
    sigcon = out_fit_no_noise{idata}.sigcon; % should be 0
    thrcon = out_fit_no_noise{idata}.thrcon;
    % simulate model
    s = stim+sigsen*randn(length(stim),nsmp);                   % stimulus sensory response according to sensory noise
    r = sign(s);                                                % model response (+1/-1)
    c = sign(abs(s)+sigcon*randn(length(stim),nsmp)-thrcon);    % model confidence (+1/-1)
    
    stim_mat = repmat(stim,1,nsmp);
    paccu_conf_tmp = [mean(r(c<0) == stim_mat(c<0));...
                      mean(r(c>0) == stim_mat(c>0))];  
    pdiff_mod_avg(idata,2) = diff(paccu_conf_tmp);
    paccu_conf_tmp = [sum((r==stim_mat)&c>0)./sum(c>0);...
                          sum((r==stim_mat)&c<0)./sum(c<0)];
    pdiff_mod_std(idata,2) = sqrt(sum(var(paccu_conf_tmp,0,2)));
    
    % per session
    for isess = 1:2
        resp = out_fit_session{idata,isess}.cfg.resp;
        stim = out_fit_session{idata,isess}.cfg.stim;
        conf = out_fit_session{idata,isess}.cfg.conf;
        pconf_sess(idata,isess) = mean(conf > 0);
        paccu_conf_tmp = [mean(resp(conf < 0) == stim(conf < 0));...  % unconfident
            mean(resp(conf > 0) == stim(conf > 0))];    % confident
        pdiff_sess_avg(idata,isess) = diff(paccu_conf_tmp);
        pdiff_sess_std(idata,isess) = sqrt( ...
            paccu_conf_tmp(1)*(1-paccu_conf_tmp(1))/nnz(conf < 0)+ ... % estimated s.d. depending on number of trials
            paccu_conf_tmp(2)*(1-paccu_conf_tmp(2))/nnz(conf > 0));
    end
    fprintf('done.\n')
    
end

% define colors
rgb = [141 191 86]/255;     % green
rgb_no_noise = [.5 .5 .5];  % grey (model with no confidence noise)

% confidence
xavg = pconf_dat; % human data
yavg = pconf_mod; % best-fitting models
% run robust regression model between human data and (full) best-fitting
% confidence model
fprintf('<strong>Robust regression model for confidence:</strong>\n');
modfun = 'y ~ b1 + b2*x1';
b0 = zeros(1,2);
opt = statset('nlinfit');
opt.RobustWgtFun = 'bisquare';
nlm = fitnlm(xavg,yavg(:,1),modfun,b0,'Options',opt);
disp(nlm);
fprintf('\n\n');
% run classic correlation
[rho,pval] = corr(xavg,yavg(:,1));
fprintf('<strong>Normal correlation for confidence:</strong>\n');
fprintf('\tRHO = %.3f ; PVAL = %.2e ; RHO^2 = %.3f\n\n',rho,pval,rho^2);
pbar = 1; HGT = 4;
figname = 'confidence_corr_dat_model_conf';
figure('Color','white','Name',figname);
hold on
xlim([0 1]);
ylim([0 1]);
plot(xlim,ylim,'k-');
plot(xavg,yavg(:,1),'wo','MarkerSize',4,'MarkerFaceColor',rgb,'LineWidth',0.5);
plot(xavg,yavg(:,2),'o','MarkerSize',3.5,'MarkerFaceColor',rgb_no_noise,'MarkerEdgeColor','non');
text(.1*diff(xlim)+min(xlim),.9*diff(ylim)+min(ylim),sprintf('r2 =  %5.3f ; p < %.2e',rho^2,pval));
set(gca,'XTick',0:.5:1,'YTick',0:.5:1);
xlabel('human data'); ylabel('best-fitting model');
helper_ILL(gcf,pbar,HGT,false);

% accuracy difference
xavg = pdiff_dat_avg; % human data
yavg = pdiff_mod_avg; % best-fitting models
xstd = pdiff_dat_std;
ystd = pdiff_mod_std;
% run robust regression model between human data and (full) best-fitting
% confidence model
fprintf('<strong>Robust regression model for accuracy difference:</strong>\n');
modfun = 'y ~ b1 + b2*x1';
b0 = zeros(1,2);
opt = statset('nlinfit');
opt.RobustWgtFun = 'bisquare';
nlm = fitnlm(xavg,yavg(:,1),modfun,b0,'Options',opt);
disp(nlm);
fprintf('\n\n');
% run classic correlation
[rho,pval] = corr(xavg,yavg(:,1));
fprintf('<strong>Normal correlation for accuracy difference:</strong>\n');
fprintf('\tRHO = %.3f ; PVAL = %.2e ; RHO^2 = %.3f\n\n',rho,pval,rho^2);
pbar = 1; HGT = 4;
figname = 'confidence_corr_dat_model_accdiff';
figure('Color','white','Name',figname);
hold on
xlim([-0.05,+0.3]);
ylim([-0.05,+0.3]);
plot(xlim,[0,0],'-','Color',[0.8,0.8,0.8]);
plot([0,0],ylim,'-','Color',[0.8,0.8,0.8]);
plot(xlim,ylim,'k-');
for i = 1:ndata
    plot(xavg(i,1)+xstd(i,1)*[+1,-1],yavg(i,1)*[1,1],'-','Color',0.5*(rgb_no_noise+1));
    plot(xavg(i,1)*[1,1],yavg(i,1)+ystd(i,1)*[+1,-1],'-','Color',0.5*(rgb_no_noise+1));
    plot(xavg(i,1)+xstd(i,1)*[+1,-1],yavg(i,2)*[1,1],'-','Color',0.5*(rgb+1));
    plot(xavg(i,1)*[1,1],yavg(i,2)+ystd(i,2)*[+1,-1],'-','Color',0.5*(rgb+1));
end
plot(xavg(:,1),yavg(:,1),'wo','MarkerSize',4,'MarkerFaceColor',rgb,'LineWidth',0.5);
plot(xavg(:,1),yavg(:,2),'o','MarkerSize',3.5,'MarkerFaceColor',rgb_no_noise,'MarkerEdgeColor','non');
text(.1*diff(xlim)+min(xlim),.9*diff(ylim)+min(ylim),sprintf('r2 =  %5.3f ; p < %.2e',rho^2,pval));
set(gca,'XTick',0:0.1:0.3,'YTick',0:0.1:0.3);
xlabel('human data'); ylabel('best-fiting model');
helper_ILL(gcf,pbar,HGT,false);

% between-session reliability : confidence
X = pconf_sess(:,1);
Y = pconf_sess(:,2);
% run robust regression model
fprintf('<strong>Robust regression model for between-session confidence:</strong>\n');
fprintf('<strong>adjusted R-Squared value reported on plot</strong>\n')
modfun = 'y ~ b1 + b2*x1';
b0 = zeros(1,2);
% opt = statset('nlinfit');
% opt.RobustWgtFun = 'bisquare';
nlm = fitnlm(X,Y,modfun,b0);%,'Options',opt);
disp(nlm);
fprintf('\n\n');
cfg = [];
cfg.rgb = rgb; % green
cfg.figname = 'confidence_between_session_conf';
cfg.xavg = X;
cfg.yavg = Y;
cfg.xlabel = 'session 1';
cfg.ylabel = 'session 2';
cfg.verbose = false;
cfg.xlim = [0 1]; cfg.xtick = 0:.5:1;
cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick;

helper_plot_corr(cfg)

% between-session reliability : accuracy difference
X = pdiff_sess_avg(:,1);
Y = pdiff_sess_avg(:,2);
X_std = pdiff_sess_std(:,1);
Y_std = pdiff_sess_std(:,2);
% run robust regression model
fprintf('<strong>Robust regression model for between-session accuracy difference:</strong>\n');
fprintf('<strong>adjusted R-Squared value reported on plot</strong>\n')
modfun = 'y ~ b1 + b2*x1';
b0 = zeros(1,2);
opt = statset('nlinfit');
opt.RobustWgtFun = 'bisquare';
nlm = fitnlm(X,Y,modfun,b0,'Options',opt);
disp(nlm);
fprintf('\n\n');
cfg = [];
cfg.rgb = rgb; % green
cfg.figname = 'confidence_between_session_accdiff';
cfg.xavg = X; cfg.xstd = X_std;
cfg.yavg = Y; cfg.ystd = Y_std;
cfg.xlabel = 'session 1';
cfg.ylabel = 'session 2';
cfg.verbose = false;
cfg.xlim = [-.05 .25]; cfg.xtick = 0:.1:0.2;
cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick;

helper_plot_corr(cfg)

% between-session reliability : best-fitting parameters
xavg_fit = cellfun(@(x)getfield(x,'xavg'),out_fit_session,'UniformOutput',0);
xstd_fit = cellfun(@(x)getfield(x,'xstd'),out_fit_session,'UniformOutput',0);
parnam = out_fit{1}.xnam;
npar   = numel(parnam);

for ipar = 1:npar
    if strcmp(parnam{ipar},'sigsen')
        continue
    end
    fprintf('Parameter %s:\n',parnam{ipar});
    
    X = cellfun(@(c)c(ipar),xavg_fit(:,1));
    Y = cellfun(@(c)c(ipar),xavg_fit(:,2));
    X_std = cellfun(@(c)c(ipar),xstd_fit(:,1));
    Y_std = cellfun(@(c)c(ipar),xstd_fit(:,2));
    % run robust regression model
    fprintf('<strong>Robust regression model for between-session best-fitting parameter %s:</strong>\n',parnam{ipar});
    fprintf('<strong>adjusted R-Squared value reported on plot</strong>\n')
    modfun = 'y ~ b1 + b2*x1';
    b0 = zeros(1,2);
    opt = statset('nlinfit');
    opt.RobustWgtFun = 'bisquare';
    nlm = fitnlm(X,Y,modfun,b0,'Options',opt);
    disp(nlm);
    fprintf('\n\n');
    cfg = [];
    cfg.rgb = rgb; % green
    cfg.figname =  sprintf('confidence_between_session_par_%s',parnam{ipar});
    cfg.xavg = X; cfg.xstd = X_std;
    cfg.yavg = Y; cfg.ystd = Y_std;
    cfg.xlabel = 'session 1';
    cfg.ylabel = 'session 2';
    cfg.verbose = false;
    
    if strcmp(parnam{ipar},'thrcon')
        cfg.xlim = [-3 5]; cfg.xtick = -2:2:4; cfg.xticklabel = {'-2.0' '0' '2.0' '4.0'};
        cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;
    elseif strcmp(parnam{ipar},'sigcon')
        cfg.xlim = [0,7.5]; cfg.xtick = 0:3:7.5; cfg.xticklabel = {'0' '3.0' '6.0'};
        cfg.ylim = cfg.xlim; cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;
    end
    
    helper_plot_corr(cfg)
    fprintf('<strong>Normal correlation for between-session best-fitting parameter %s:</strong>\n',parnam{ipar});
    [rho,pval]=corr(cfg.xavg ,cfg.yavg);
    [R,P,RLO,RUP]=corrcoef(cfg.xavg ,cfg.yavg);
    fprintf('\tRHO = %.3f ; PVAL = %.3f ; RHO^2 = %.3f ; CI = [%+.3f %+.3f]\n\n',rho,pval,rho^2,RLO(1,2),RUP(1,2));
end

% Reliability of positive relation between reliability threshold and confidence threshold
% set color for conditional inference model
rgb = [134,198,226]/255;
% get best-fitting confidence threshold parameter:
load(fullfile(savepath,'expe2_results_fit_confidence.mat'),'out_fit');
xavg_fit = cell2mat(cellfun(@(s) s.xavg ,out_fit,'UniformOutput',false));
xstd_fit = cell2mat(cellfun(@(s) s.xstd ,out_fit,'UniformOutput',false));
sigcon_avg = xavg_fit(:,2);
sigcon_std = xstd_fit(:,2);
thrcon_avg = xavg_fit(:,3);
thrcon_std = xstd_fit(:,3);
% get best-fitting reliability threshold parameter:
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');
imodel = 3; % conditional inference model
xavg_fit = cell2mat(cellfun(@(s) s.xavg ,out_fit(:,imodel),'UniformOutput',false));
xstd_fit = cell2mat(cellfun(@(s) s.xstd ,out_fit(:,imodel),'UniformOutput',false));
delta_avg = xavg_fit(:,3);
delta_std = xstd_fit(:,3);

% run multiple regression model between confidence parameters and
% conditional inference threshold
[~,~,stat] = glmfit(zscore([thrcon_avg sigcon_avg],[],1),delta_avg,'normal'); % z-score to have comparable reg coefficients
fprintf('<strong>Multiple regression model:</strong>\n');
fprintf('  * b(thrcon) = %+.3f +/- %.3f ; p = %.3f ; t = %.3f ; CI = [%+.3f %+.3f]\n',stat.beta(2),stat.se(2),stat.p(2),stat.t(2),stat.beta(2)-1.96*stat.se(2),stat.beta(2)+1.96*stat.se(2));
fprintf('  * b(sigcon) = %+.3f +/- %.3f ; p = %.3f ; t = %.3f ; CI = [%+.3f %+.3f]\n',stat.beta(3),stat.se(3),stat.p(3),stat.t(3),stat.beta(3)-1.96*stat.se(3),stat.beta(3)+1.96*stat.se(3));
fprintf('\n');


pbar = 1;
% define predictions and confidence intervals
for ipar = 1:2
    [b,~,stat] = glmfit([thrcon_avg sigcon_avg],delta_avg,'normal');
    if ipar == 1 % thrcon
        figure('Color','white','Name','confidence_multiple_regression_thrcon')
        X = thrcon_avg;
        X_std = thrcon_std;
        xvs = linspace(-3,+5,100);
        [yvs,yci(1,:),yci(2,:)] = glmval(b, ...
            [xvs(:) mean(sigcon_avg)*ones(100,1)],'identity',stat);
        xt = -2:2:4; xtl = {'-2.0' '0' '2.0' '4.0'};xl = 'confidence threshold';
    else % sigcon
        figure('Color','white','Name','confidence_multiple_regression_sigcon')
        X = sigcon_avg;
        X_std = sigcon_std;
        xvs = linspace(0,7.5,100);
        [yvs,yci(1,:),yci(2,:)] = glmval(b, ...
            [mean(thrcon_avg)*ones(100,1),xvs(:)],'identity',stat);
        xt = 0:3:6; xtl = {'0' '3.0' '6.0'}; xl = 'confidence noise';
    end
    yci = [yvs'-yci(1,:);yvs'+yci(2,:)];
    Y = delta_avg; Y_std = delta_std;
    hold on
    patch([xvs,fliplr(xvs)],[yci(1,:),fliplr(yci(2,:))],0.5*(rgb+1),'EdgeColor','none','FaceAlpha',0.5);
    plot(xvs,yvs,'-','Color',rgb,'LineWidth',2);
    for i = 1:ndata
        plot(X(i)+X_std(i)*[+1,-1],Y(i)*[1,1],'-','Color',0.5*(rgb+1));
        plot(X(i)*[1,1],Y(i)+Y_std(i)*[+1,-1],'-','Color',0.5*(rgb+1));
    end
    plot(X,Y,'wo','MarkerSize',4,'MarkerFaceColor',rgb,'LineWidth',0.5);
    ylim([0,2.5]); xlim([min(xvs) max(xvs)]);
    set(gca,'XTick',xt,'XTickLabel',xtl,'YTick',0:1:2.5,'YTickLabel',{'0' '1.0' '2.0'});
    ylabel('reliability threshold');xlabel(xl)
    helper_ILL(gcf,pbar,HGT,false);
end

%% SUPPLEMENTARY: BETWEEN-SESSION RELIABILITY (Supplementary Figure 6)
fprintf('<strong>*** SUPPLEMENTARY: BETWEEN-SESSION RELIABILITY ***</strong>\n\n');

clearvars;
savepath = '../results';

rgb = [ ...
    238,197,141; ... % senbias
    048,102,158; ... % inflaps
    134,198,226; ... % infdisc
    171,131,171; ... % selbias
    099,074,135  ... % selepsi
    ]/255;

load(fullfile(savepath,'expe1_results_fit_stabilizing_models_per_session.mat'),'out_fit')
out_fit_expe1 = out_fit;
ndata1 = size(out_fit_expe1,1);
load(fullfile(savepath,'expe2_results_fit_stabilizing_models_per_session.mat'),'out_fit')
out_fit = cat(1,out_fit_expe1,out_fit);

% conditional inference model
imodel = 3; % conditional inference
% get best-fitting parameters
xavg_fit = cell2mat(cellfun(@(x)getfield(x,'xavg'),out_fit(:,imodel,:),'UniformOutput',0));
xstd_fit = cell2mat(cellfun(@(x)getfield(x,'xstd'),out_fit(:,imodel,:),'UniformOutput',0));
parnam = out_fit{1,imodel,1}.xnam;
npar   = numel(parnam);

for ipar = 1:npar
    fprintf('Parameter %s:\n',parnam{ipar});
    
    X = xavg_fit(:,ipar,1);
    Y = xavg_fit(:,ipar,2);
    X_std = xstd_fit(:,ipar,1);
    Y_std = xstd_fit(:,ipar,2);
    
    cfg = [];
    cfg.rgb = rgb(imodel,:);
    cfg.figname =  sprintf('between_session_par_%s',parnam{ipar});
    cfg.xavg = X; cfg.xstd = X_std;
    cfg.yavg = Y; cfg.ystd = Y_std;
    cfg.xlabel = 'session 1';
    cfg.ylabel = 'session 2';
    
    if strcmp(parnam{ipar},'h')
        cfg.xlim = [0 .4]; cfg.xtick = 0:.1:.4; cfg.xticklabel = {'0' '0.1' '0.2' '0.3' '0.4'};
        cfg.ylim = [0 .4]; cfg.ytick = 0:.1:.4; cfg.yticklabel = {'0' '0.1' '0.2' '0.3' '0.4'};
    elseif strcmp(parnam{ipar},'siginf')
        cfg.xlim = [0 4]; cfg.xtick = 0:1:4; cfg.xticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
        cfg.ylim = [0 4]; cfg.ytick = 0:1:4; cfg.yticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(parnam{ipar},'delta')
        cfg.xlim = [0 3.5]; cfg.xtick = 0:1.0:3.5; cfg.xticklabel = {'0' '1.0' '2.0' '3.0'};
        cfg.ylim = [0 3.5]; cfg.ytick = 0:1.0:3.5; cfg.yticklabel = {'0' '1.0' '2.0' '3.0'};
    elseif strcmp(parnam{ipar},'beta')
        cfg.xlim = [0 1.5]; cfg.xtick = 0:.5:1.5; cfg.xticklabel = {'0' '0.5' '1.0' '1.5'};
        cfg.ylim = [0 1.5]; cfg.ytick = 0:.5:1.5; cfg.yticklabel = {'0' '0.5' '1.0' '1.5'};
    elseif strcmp(parnam{ipar},'epsi')
        cfg.xlim = [0 .1]; cfg.xtick = 0:.05:0.1; cfg.xticklabel = {'0' '0.05' '0.1'};
        cfg.ylim = [0 .1]; cfg.ytick = 0:.05:0.1; cfg.yticklabel = {'0' '0.05' '0.1'};
    end
    
    helper_plot_corr(cfg)
    [~,~,RLO,RUP]=corrcoef(cfg.xavg ,cfg.yavg);
    fprintf('\t\tCI = [%+.3f %+.3f]\n\n',RLO(1,2),RUP(1,2));
end

% Between-session effect on response switch and reversal curves
crev_dat = cat(1,squeeze(cellfun(@(s)getfield(s,'crev_sub'),out_fit(1:ndata1,imodel,:),'UniformOutput',0)),...
           cellfun(@(x)mean(x,2),squeeze(cellfun(@(x)getfield(x,'crev_w_sub'),out_fit(ndata1+1:end,imodel,:),'UniformOutput',0)),'UniformOutput',0));
crep_dat = cat(1,squeeze(cellfun(@(s)getfield(s,'crep_sub'),out_fit(1:ndata1,imodel,:),'UniformOutput',0)),...
           cellfun(@(x)mean(x,2),squeeze(cellfun(@(x)getfield(x,'crep_w_sub'),out_fit(ndata1+1:end,imodel,:),'UniformOutput',0)),'UniformOutput',0));
crev_mod = cat(1,squeeze(cellfun(@(s)getfield(s,'crev_avg'),out_fit(1:ndata1,imodel,:),'UniformOutput',0)),...
           cellfun(@(x)mean(x,2),squeeze(cellfun(@(x)getfield(x,'crev_w_avg'),out_fit(ndata1+1:end,imodel,:),'UniformOutput',0)),'UniformOutput',0));
crep_mod = cat(1,squeeze(cellfun(@(s)getfield(s,'crep_avg'),out_fit(1:ndata1,imodel,:),'UniformOutput',0)),...
           cellfun(@(x)mean(x,2),squeeze(cellfun(@(x)getfield(x,'crep_w_avg'),out_fit(ndata1+1:end,imodel,:),'UniformOutput',0)),'UniformOutput',0));

crev_dat = cellfun(@(x) x',crev_dat,'UniformOutput',0);
crep_dat = cellfun(@(x) x',crep_dat,'UniformOutput',0);
crev_mod = cellfun(@(x) x',crev_mod,'UniformOutput',0);
crep_mod = cellfun(@(x) x',crep_mod,'UniformOutput',0); 

serr = @(x,d)std(x,[],d)/sqrt(size(x,d)); % standard error

for ipar = 1:npar
    % session-wise
    % mediansplit on session 1
    X = xavg_fit(:,ipar,1);
    Xmed = median(X);
    idx1 = X<=Xmed;
    % mediansplit on session 2
    X = xavg_fit(:,ipar,2);
    Xmed = median(X);
    idx2 = X<=Xmed;
    
    % mean+serr crossvalidated REV+REP curves
    crev_dat_avg(1,:) = mean(cell2mat(cat(1,crev_dat(idx1,2),crev_dat(idx2,1))));
    crev_dat_avg(2,:) = mean(cell2mat(cat(1,crev_dat(~idx1,2),crev_dat(~idx2,1))));
    crep_dat_avg(1,:) = mean(cell2mat(cat(1,crep_dat(idx1,2),crep_dat(idx2,1))));
    crep_dat_avg(2,:) = mean(cell2mat(cat(1,crep_dat(~idx1,2),crep_dat(~idx2,1))));
    crev_mod_avg(1,:) = mean(cell2mat(cat(1,crev_mod(idx1,2),crev_mod(idx2,1))));
    crev_mod_avg(2,:) = mean(cell2mat(cat(1,crev_mod(~idx1,2),crev_mod(~idx2,1))));
    crep_mod_avg(1,:) = mean(cell2mat(cat(1,crep_mod(idx1,2),crep_mod(idx2,1))));
    crep_mod_avg(2,:) = mean(cell2mat(cat(1,crep_mod(~idx1,2),crep_mod(~idx2,1))));
    
    crev_dat_err(1,:) = serr(cell2mat(cat(1,crev_dat(idx1,2),crev_dat(idx2,1))),1);
    crev_dat_err(2,:) = serr(cell2mat(cat(1,crev_dat(~idx1,2),crev_dat(~idx2,1))),1);
    crep_dat_err(1,:) = serr(cell2mat(cat(1,crep_dat(idx1,2),crep_dat(idx2,1))),1);
    crep_dat_err(2,:) = serr(cell2mat(cat(1,crep_dat(~idx1,2),crep_dat(~idx2,1))),1);
    crev_mod_err(1,:) = serr(cell2mat(cat(1,crev_mod(idx1,2),crev_mod(idx2,1))),1);
    crev_mod_err(2,:) = serr(cell2mat(cat(1,crev_mod(~idx1,2),crev_mod(~idx2,1))),1);
    crep_mod_err(1,:) = serr(cell2mat(cat(1,crep_mod(idx1,2),crep_mod(idx2,1))),1);
    crep_mod_err(2,:) = serr(cell2mat(cat(1,crep_mod(~idx1,2),crep_mod(~idx2,1))),1);
    
    % SWITCH curve
    cfg = [];
    cfg.figname = sprintf('between_session_SWI_%s',parnam{ipar});
    cfg.rgb = [(1+rgb(imodel,:))/2;rgb(imodel,:)];
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-crep_mod_avg';
    cfg.yerr = crep_mod_err';
    cfg.Yavg = 1-crep_dat_avg';
    cfg.Yerr = crep_dat_err';
    helper_plot_rep(cfg)
    
    % REVERSAL curve
    cfg = [];
    cfg.figname = sprintf('between_session_REV_%s',parnam{ipar});
    cfg.rgb = [(1+rgb(imodel,:))/2;rgb(imodel,:)];
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = crev_mod_avg';
    cfg.yerr = crev_mod_err';
    cfg.Yavg = crev_dat_avg';
    cfg.Yerr = crev_dat_err';
    helper_plot_rev(cfg)
    
end

% variance explained
vexpl_dat = nan(3,2,2); % variance explained for reversal=1 and repetition=2 curves for human data
modelfun = @(b,x)b(1)+b(2)*x; % model function
start = [0;0]; % starting values
fprintf('\n<strong>Variance explained by conditional inference model</strong>\n');
for ipar = 1:npar
    X = xavg_fit(:,ipar,1); % fitted parameters SESSION 1
    xdat = cat(2,cell2mat(crev_dat(:,2)),cell2mat(crep_dat(:,2))); % subject data SESSION 2
    xres = nan(size(xdat));
    r2 = nan(1,16);
    for j = 1:16
        Y = xdat(:,j);
        nlm = fitnlm(X,Y,modelfun,start);
        xres(:,j) = nlm.Residuals.Raw;
        r2(j) = corr(X,Y)^2;
    end
    vexpl_dat(ipar,1,1) = 1-sum(var(xres(:,1:8)))/sum(var(xdat(:,1:8)));    % REV
    vexpl_dat(ipar,2,1) = 1-sum(var(xres(:,9:16)))/sum(var(xdat(:,9:16)));  % REP

    X = xavg_fit(:,ipar,2); % fitted parameters SESSION 2
    xdat = cat(2,cell2mat(crev_dat(:,1)),cell2mat(crep_dat(:,1))); % subject data SESSION 1
    xres = nan(size(xdat));
    for j = 1:16
        Y = xdat(:,j);
        nlm = fitnlm(X,Y,modelfun,start);
        xres(:,j) = nlm.Residuals.Raw;
    end
    vexpl_dat(ipar,1,2) = 1-sum(var(xres(:,1:8)))/sum(var(xdat(:,1:8)));    % REV
    vexpl_dat(ipar,2,2) = 1-sum(var(xres(:,9:16)))/sum(var(xdat(:,9:16)));  % REP
    vexpl_dat = mean(vexpl_dat,3);
    fprintf('Parameter %s:\n',parnam{ipar})
    fprintf('\t* switch curve:   %3.1f%%\n',vexpl_dat(ipar,2)*100);
    fprintf('\t* reversal curve: %3.1f%%\n\n',vexpl_dat(ipar,1)*100);
end

% Bayesian Model Selection
elbo_session = cellfun(@(s)getfield(s,'elbo'),out_fit);

for isession = 1:2
    cfg = [];
    cfg.figname = sprintf('BMS_session%d',isession);
    cfg.rgb = rgb;
    cfg.elbo = elbo_session(:,:,isession);
    cfg.ylabel = 'model probability';
    cfg.xlabel = 'model';
    helper_plot_bms(cfg);
end

%% SUPPLEMENTARY: PRINCIPAL COMPONENT ANALYSIS (Supplementary Figure7)
fprintf('<strong>*** SUPPLEMENTARY: PRINCIPAL COMPONENT ANALYSIS ***</strong>\n\n');

clearvars;
savepath = '../results';

% define useful function handles
serr = @(x,d)std(x,[],d)/sqrt(size(x,d)); % s.e.m.
shuffle = @(x)x(randperm(numel(x))); % shuffle vector

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapse
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapse
    ]/255;

% load stabilizing models fitting results from both experiments
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit')
out_fit_expe1 = out_fit;
ndata1 = size(out_fit_expe1,1);
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit')
out_fit = cat(1,out_fit_expe1,out_fit);

% conditional inference model
imodel = 3; % conditional inference
% get best-fitting parameters
xavg_fit = cell2mat(cellfun(@(x)getfield(x,'xavg'),out_fit(:,imodel),'UniformOutput',0));
xstd_fit = cell2mat(cellfun(@(x)getfield(x,'xstd'),out_fit(:,imodel),'UniformOutput',0));
parnam = out_fit{1,imodel}.xnam;
npar   = numel(parnam);
% get reversal metrics (human data and conditional inference model)
crev_dat = cat(1,cellfun(@(s)getfield(s,'crev_sub'),out_fit(1:ndata1,imodel),'UniformOutput',0),...
           cellfun(@(x)mean(x,2),cellfun(@(x)getfield(x,'crev_w_sub'),out_fit(ndata1+1:end,imodel),'UniformOutput',0),'UniformOutput',0));
crep_dat = cat(1,cellfun(@(s)getfield(s,'crep_sub'),out_fit(1:ndata1,imodel),'UniformOutput',0),...
           cellfun(@(x)mean(x,2),cellfun(@(x)getfield(x,'crep_w_sub'),out_fit(ndata1+1:end,imodel),'UniformOutput',0),'UniformOutput',0));
crev_mod = cat(1,cellfun(@(s)getfield(s,'crev_avg'),out_fit(1:ndata1,imodel),'UniformOutput',0),...
           cellfun(@(x)mean(x,2),cellfun(@(x)getfield(x,'crev_w_avg'),out_fit(ndata1+1:end,imodel),'UniformOutput',0),'UniformOutput',0));
crep_mod = cat(1,cellfun(@(s)getfield(s,'crep_avg'),out_fit(1:ndata1,imodel,:),'UniformOutput',0),...
           cellfun(@(x)mean(x,2),cellfun(@(x)getfield(x,'crep_w_avg'),out_fit(ndata1+1:end,imodel),'UniformOutput',0),'UniformOutput',0));
crev_dat = cell2mat(cellfun(@(x) x',crev_dat,'UniformOutput',0));
crep_dat = cell2mat(cellfun(@(x) x',crep_dat,'UniformOutput',0));
crev_mod = cell2mat(cellfun(@(x) x',crev_mod,'UniformOutput',0));
crep_mod = cell2mat(cellfun(@(x) x',crep_mod,'UniformOutput',0));

ndata = size(crev_dat,1);

% set data arrays
xdat = cat(2,crev_dat,crep_dat); % human data
xmod = cat(2,crev_mod,crep_mod); % best-fitting conditional inference model predictions

% compute SVD-based PCA on human data
[~,score_pca,latnt_pca] = pca(xdat);

npc = npar;         % number of principal components to consider
npt = size(xdat,2); % number of data points
% get variance explained per principal component
vexpl_tot = latnt_pca/sum(latnt_pca); % total variance explained
vexpl_dat = nan(npc,2); % variance explained for reversal=1 and repetition=2 curves for human data
vexpl_mod = nan(npc,2); % variance explained for reversal=1 and repetition=2 curves for model predictions

for ipc = 1:npc
    fprintf('<strong>Variance explained by PC%d</strong>\n',ipc)
    % participants
    xres = nan(ndata,npt);
    for j = 1:npt
        b = [ones(ndata,1),score_pca(:,ipc)]\xdat(:,j);
        xres(:,j) = xdat(:,j)-[ones(ndata,1),score_pca(:,ipc)]*b;
    end
    vexpl_dat(ipc,1) = 1-sum(var(xres(:,1:8)))/sum(var(xdat(:,1:8)));
    vexpl_dat(ipc,2) = 1-sum(var(xres(:,9:16)))/sum(var(xdat(:,9:16)));
    % model
    xres = nan(ndata,npt);
    for j = 1:npt
        b = [ones(ndata,1),score_pca(:,ipc)]\xmod(:,j);
        xres(:,j) = xmod(:,j)-[ones(ndata,1),score_pca(:,ipc)]*b;
    end
    vexpl_mod(ipc,1) = 1-sum(var(xres(:,1:8)))/sum(var(xmod(:,1:8)));
    vexpl_mod(ipc,2) = 1-sum(var(xres(:,9:16)))/sum(var(xmod(:,9:16)));
    fprintf('\tTOT: %02.1f%%\n\tSWI: %02.1f%%\n\tREV: %02.1f%%\n\n',vexpl_tot(ipc)*100,vexpl_dat(ipc,2)*100,vexpl_dat(ipc,1)*100)
    
    % plot switch and reversal curves
    [rho,pval] = corr(score_pca(:,ipc),xavg_fit(:,ipc));
    X = score_pca(:,ipc);
    if rho < 0
        X = -score_pca(:,ipc);
    end
    xmed = median(X);
    idx = X <= xmed;
    
    % RESPONSE SWITCH curves
    crep_mod_avg = cat(2,mean(crep_mod(idx,:))',mean(crep_mod(~idx,:))');
    crep_mod_err = cat(2,serr(crep_mod(idx,:),1)',serr(crep_mod(~idx,:),1)');
    crep_dat_avg = cat(2,mean(crep_dat(idx,:))',mean(crep_dat(~idx,:))');
    crep_dat_err = cat(2,serr(crep_dat(idx,:),1)',serr(crep_dat(~idx,:),1)');
    cfg = [];
    cfg.figname = sprintf('PCA_SWI_PC%d',ipc);
    cfg.rgb = [(1+rgb(imodel,:))/2;rgb(imodel,:)];
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-crep_mod_avg;
    cfg.yerr = crep_mod_err;
    cfg.Yavg = 1-crep_dat_avg;
    cfg.Yerr = crep_dat_err;
    helper_plot_rep(cfg)
    
    % RESPONSE REVERSAL curves
    crev_mod_avg = cat(2,mean(crev_mod(idx,:))',mean(crev_mod(~idx,:))');
    crev_mod_err = cat(2,serr(crev_mod(idx,:),1)',serr(crev_mod(~idx,:),1)');
    crev_dat_avg = cat(2,mean(crev_dat(idx,:))',mean(crev_dat(~idx,:))');
    crev_dat_err = cat(2,serr(crev_dat(idx,:),1)',serr(crev_dat(~idx,:),1)');
    cfg = [];
    cfg.figname = sprintf('PCA_REV_PC%d',ipc);
    cfg.rgb = [(1+rgb(imodel,:))/2;rgb(imodel,:)];
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = crev_mod_avg;
    cfg.yerr = crev_mod_err;
    cfg.Yavg = crev_dat_avg;
    cfg.Yerr = crev_dat_err;
    helper_plot_rev(cfg)
    
end

% get correlation between each PC score and each model parameter
nres   = 1e4;               % number of bootstrap resamples
rsq    = nan(npc,npar);     % r-squared
rsq_se = nan(npc,npar);     % bootstrap estimate of standard error
rsq_qt = nan(npc,npar,2);   % bootstrap estimates of 1st/3rd quartiles
rsq_ci = nan(npc,npar,2);   % bootstrap estimates of 95% CI
pexc   = nan(npc,npar);     % exceedance probab
rho    = nan(npc,npar);
for ipc = 1:npc
    for j = 1:npar
        [rho(ipc,j),pval] = corr(score_pca(:,ipc),xavg_fit(:,j));
        rsq(ipc,j) = rho(ipc,j)^2;
        if ipc == j
            fprintf('<strong>PC%d score versus %s</strong>\n',ipc,parnam{j})
            fprintf('\tr = %+.3f, p = %.3f, r^2 = %+.3f\n',abs(rho(ipc,j)),pval,rho(ipc,j)^2);
        end
    end
    rsq_res = nan(nres,3);
    for ires = 1:nres
        isub = randsample(ndata,ndata,true);
        for j = 1:npar
            rho_tmp = corr(score_pca(isub,ipc),xavg_fit(isub,j));
            rsq_res(ires,j) = rho_tmp^2;
        end
    end
    [~,jmax] = max(rsq_res,[],2);
    for j = 1:npar
        pexc(ipc,j) = mean(jmax == j);
        rsq_se(ipc,j) = std(rsq_res(:,j));
        rsq_qt(ipc,j,:) = quantile(rsq_res(:,j),[1,3]/4);
        rsq_ci(ipc,j,:) = quantile(rsq_res(:,j),[0.025 0.975]);
        if ipc == j
            fprintf('\texceedance p = %.3f, 95%% CI = [%+.3f %+.3f]\n\n',pexc(ipc,j),rsq_ci(ipc,j,1),rsq_ci(ipc,j,2));
        end
    end
    
    % plot variance of each principal component explained by corresponding model parameter
    cfg = [];
    cfg.figname = sprintf('PC%d_variance_explained',ipc);
    [~,Imax] = max(pexc(ipc,:));
    cfg.rgb = repmat((1+rgb(imodel,:))/2,npar,1);
    cfg.rgb(Imax,:) = rgb(imodel,:);
    cfg.x = rsq(ipc,:);
    cfg.xqt =  squeeze(rsq_qt(ipc,:,:));
    if max(squeeze(rsq_qt(ipc,:,:)),[],'all')<0.5
        cfg.ylim = [0 .5];
    else
        cfg.ylim = [0 1];
    end
    cfg.ylabel = 'variance explained';
    cfg.xlabel = 'model parameter';
    cfg.xticklabel = parnam;
    cfg.HGT = 4; cfg.pbar = 1;
    helper_plot_overall_bars(cfg);
    
    % plot correlation between conditional inference model parameter and principal component
    if rho(ipc,ipc)<0
        X = -score_pca(:,ipc);
    else
        X = +score_pca(:,ipc);
    end
    cfg = [];
    cfg.rgb = rgb(imodel,:);
    cfg.xavg = X;
    cfg.yavg = xavg_fit(:,ipc);
    cfg.ystd = xstd_fit(:,ipc);
    cfg.ci = true;
    if ipc == 1
        xlim = [-.8,.6]; ylim = [0 .4];
        xtick = -.8:.4:.6; ytick = 0:.1:.4;
    elseif ipc == 2
        xlim = [-.6,.4]; ylim = [0 4];
        xtick = -4:.4:.4; ytick = 0:1:4;
    elseif ipc == 3
        xlim = [-.4,.4]; ylim = [0 3];
        xtick = -.4:.2:.4; ytick = 0:1:3;
    end
    cfg.xlim = xlim; cfg.ylim = ylim;
    cfg.xtick = xtick; cfg.ytick = ytick;
    cfg.xlabel = sprintf('PC%d score',ipc);
    cfg.ylabel = parnam{ipc};
    cfg.verbose = false;
    helper_plot_corr(cfg);
end

%% SUPPLEMENTARY: CONDITIONAL INFERENCE MODEL EFFECTS (Supplementary Figure 8)
fprintf('<strong>*** SUPPLEMENTARY: CONDITIONAL INFERENCE MODEL EFFECTS ***</strong>\n\n');

clearvars;
savepath = '../results';
% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapses
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapses
    ]/255;

% load simulations output
load(fullfile(savepath,'expe_merged_predicted_effects_accuracy.mat'));
% plot predicted accuracy on the last trial before and the first trial
% after a reversal varying each stabilization parameter

pbar = .9; HGT = 4;
for imodel = 1:nmodel
    xvec = stab_param_vec(imodel,:);
    Y1 = squeeze(mean(firstaf_sim_avg(:,:,imodel)));
    Y2 = 1-squeeze(mean(lastbef_sim_avg(:,:,imodel)));
    pcor0_1 = interp1(xvec,Y1,0,'spline');
    pcor0_2 = interp1(xvec,Y2,0,'spline');    
    % last before / first after
    figname = sprintf('accuracy_before_after_model_%s_%s',modnam{imodel},stab_param_nam{imodel});
    figure('Color','white','Name',figname)
    hold on
    xlim([min(xvec) max(xvec)]);
    ylim([0 1]);
    plot(xvec,Y1,'Color',rgb(imodel,:),'LineWidth',1.5);
    plot(xvec,Y2,'Color',rgb(imodel,:),'LineWidth',1.5,'LineStyle',':');
    plot(xlim,pcor0_1*[1 1],'Color',[.8 .8 .8])
    plot(xlim,pcor0_2*[1 1],'Color',[.8 .8 .8])
    switch stab_param_nam{imodel}
        case 'lambda'
            set(gca,'XTick',0:.5:1.5,'XTickLabel',{'0','0.5','1.0','1.5'});
        case {'delta','beta'}
            set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        case {'plaps','epsi'}
            set(gca,'XTick',0:.5:1);
    end
    xlabel(sprintf('%s',stab_param_nam{imodel}));
    helper_ILL(gcf,pbar,HGT,false);
end 

% Psychometric effects of conditional inference
fprintf('Psychometric effects of conditional inference\n')
imodel = 3; % conditional inference model
pbar = .9; HGT = 4;
% plot predicted accuracy for conditional inference model varying delta
xdat = stab_param_vec(imodel,:)';
ydat = mean(pcor_sim_avg(:,:,imodel),1)';
par_opt = fminbnd(@(x)-interp1(xdat,ydat,x,'spline'),0,max(xdat));
pcor_opt = interp1(xdat,ydat,par_opt,'spline');
pcor0 = interp1(xdat,ydat,0,'spline');
figname = sprintf('accuracy_simulation_model_%s',modnam{imodel});
figure('Color','white','Name',figname)
hold on
xlim([min(xdat) max(xdat)]);
ylim([0.5 0.9]);
plot(par_opt*[1,1],ylim,'k-');
plot(xdat,ydat,'-','LineWidth',1.5,'Color',rgb(imodel,:));
plot(par_opt,pcor_opt,'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
gain = (pcor_opt-pcor0)*100;
fprintf('\tOptimal accuracy:\t%4.2f%% (= p0%+4.3f) for %s \t= %4.3f\n',pcor_opt*100,gain,stab_param_nam{imodel},par_opt)
set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
plot(xlim,pcor0*[1 1],'k:')
xlabel(sprintf('%s',stab_param_nam{imodel}));
set(gca,'YTick',0.5:0.1:0.9);
ylabel('fraction correct');
helper_ILL(gcf,pbar,HGT);


% plot corresponding fraction discarded
delta_vec = xdat;
sigsen = 1/norminv(1-0.20); % sensory noise
ph_vec= 1-normcdf(+delta_vec,2/sigsen^2,2/sigsen);  % fraction hits
pf_vec = normcdf(-delta_vec,2/sigsen^2,2/sigsen);   % fraction false-alarms
pd_vec = 1-ph_vec-pf_vec;                           % fraction discarded
ydat = pd_vec;
pd_opt = interp1(xdat,ydat,par_opt,'spline');
fprintf('\tFraction discarded:\t%4.2f%% for %s \t= %4.3f\n',pd_opt*100,stab_param_nam{imodel},par_opt)
figname = 'fraction_discarded_simulation_model_infdisc';
figure('Color','white','Name',figname)
hold on
xlim([min(xdat) max(xdat)]);
plot(par_opt*[1,1],ylim,'k-');
plot(xdat,ydat,'-','LineWidth',1.5,'Color',rgb(imodel,:));
plot(par_opt,pd_opt,'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
xlabel(sprintf('%s',stab_param_nam{imodel}));
set(gca,'YTick',0:0.2:1,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
ylabel('fraction discarded');
helper_ILL(gcf,pbar,HGT);

% plot corresponding fraction switch on trials where switch is correct or
% incorrect
xdat = stab_param_vec(imodel,:);
ydat = squeeze(mean(pswi_sim_avg(:,:,imodel,:),1));
pswi_opt = interp1(xdat,ydat,par_opt,'spline');
pswi0 = interp1(xdat,ydat,0,'spline');
figname = sprintf('fraction_switch_simulation_model_%s',modnam{imodel});
figure('Color','white','Name',figname)
pbar = .9; HGT = 4;
hold on
xlim([min(xdat) max(xdat)]);
ylim([0 .6]);
plot(par_opt*[1,1],ylim,'k-');
plot(par_opt,pswi_opt(1),'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
plot(par_opt,pswi_opt(2),'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
plot(xdat,ydat(:,1),'Color',rgb(imodel,:),'LineWidth',1.5);
plot(xdat,ydat(:,2),'Color',rgb(imodel,:),'LineWidth',1.5,'LineStyle',':');
plot(xlim,pswi0(1)*[1 1],'k:');plot(xlim,pswi0(2)*[1 1],'k:');
set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
xlabel(sprintf('%s',stab_param_nam{imodel}));
helper_ILL(gcf,pbar,HGT);

% Observed and predicted effects of inference lapses model parameters on accuracy
imodel = 2; % inference lapses model
nsmp = 1e4; % number of samples for estimating peak value and s.e.

mtyp = 'd';
ydat = pcor_mod_avg(:,imodel);
ystd = pcor_mod_std(:,imodel);
xdat = cell2mat(xavg_fit(:,imodel));
xstd = cell2mat(xstd_fit(:,imodel));
% fit regression model
modfun = 'y ~ b1 + b2*x1 + b3*x1^2 + b4*x2 + b5*x3 + b6*x3^2';
b0 = cat(2,mean(ydat),zeros(1,5));
opt = statset('nlinfit');
opt.RobustWgtFun = 'bisquare';
nlm = fitnlm(xdat,ydat,modfun,b0,'Options',opt);
w = nlm.Robust.Weights;
fprintf('<strong>Regression model for fraction correct inference lapses model:</strong>\n');
disp(nlm);

% estimate performance-maximizing parameter values
fprintf('\nComputing estimate of performance-maximizing parameter values ...\n')
bsmp = mvnrnd(nlm.Coefficients.Estimate',nlm.CoefficientCovariance,nsmp);
fun = @(x,b)b(1)+b(2)*x(1)+b(3)*x(1)^2 + b(4)*x(2) + b(5)*x(3)+b(6)*x(3)^2;
options = optimoptions('fmincon','Display','none');
xpeak = nan(nsmp,3);
ppeak = nan(nsmp,1);
for ismp = 1:nsmp
    xpeak(ismp,:) = fmincon(@(x)-fun(x,bsmp(ismp,:)), ...
        mean(xdat,1)',[],[],[],[],[0;0;0],[0.4;4.0;3.0],[],options);
    ppeak(ismp) = fun(xpeak(ismp,:),bsmp(ismp,:));
end
xpeak_avg = mean(xpeak,1);
xpeak_se = std(xpeak,[],1);

for ipar = 1:npar
    fprintf('\t%s\t= %5.3f +/- %5.3f\n',parnam{imodel}{ipar},xpeak_avg(ipar),xpeak_se(ipar));
    % plot fraction correct wrt parameter ipar
    switch parnam{imodel}{ipar}
        case 'h'
            x = linspace(0,0.4,100)';
            xnew = cat(2,x,mean(xdat(:,2))*ones(100,1),mean(xdat(:,3))*ones(100,1));
        case 'siginf'
            x = linspace(0,4,100)';
            xnew = cat(2,mean(xdat(:,1))*ones(100,1),x,mean(xdat(:,3))*ones(100,1));
        case stab_param_nam{imodel}
            x = linspace(min(stab_param_vec(imodel,:)),max(stab_param_vec(imodel,:)),100)';
            xnew = cat(2,mean(xdat(:,1))*ones(100,1),mean(xdat(:,2))*ones(100,1),x);
    end
    [yhat,y_ci] = predict(nlm,xnew);
    
    figname = sprintf('accuracy_multiple_regression_model_%s_%s',modnam{imodel},parnam{imodel}{ipar});
    figure('Color','white','Name',figname)
    pbar = 1.3; HGT = 4;
    hold on
    xlim([min(x) max(x)]);
    ylim([0.6 0.9]);
    
    patch([x;flipud(x)],[y_ci(:,1);flipud(y_ci(:,2))],0.5*(rgb(imodel,:)+1),'EdgeColor','none','FaceAlpha',0.5);
    plot(xpeak_avg(ipar)*[1,1],ylim,'k-');
    plot(x,yhat,'-','LineWidth',1.5,'Color',rgb(imodel,:));
    
    for idata = 1:ndata
        rgb_w = rgb(imodel,:)*w(idata)+(1-w(idata));
        plot(xdat(idata,ipar)+xstd(idata,ipar)*[+1,-1],ydat(idata)*[1,1],'-','Color',0.5*(rgb_w+1));
        plot(xdat(idata,ipar)*[1,1],ydat(idata)+ystd(idata)*[+1,-1],'-','Color',0.5*(rgb_w+1));
    end
    
    for idata = 1:ndata
        rgb_w = rgb(imodel,:)*w(idata)+(1-w(idata));
        plot(xdat(idata,ipar),ydat(idata),['w',mtyp],'MarkerSize',4,'MarkerFaceColor',rgb_w);
    end
    switch parnam{imodel}{ipar}
        case 'h'
            set(gca,'XTick',0:0.1:0.4);
        case 'siginf'
            set(gca,'XTick',0:1:4,'XTickLabel',{'0','1.0','2.0','3.0','4.0'});
        case 'lambda'
            set(gca,'XTick',0:.5:1.5,'XTickLabel',{'0','0.5','1.0','1.5'});
        case {'delta','beta'}
            set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        case {'plaps','epsi'}
            set(gca,'XTick',0:.5:1);
    end
    xlabel(sprintf('%s',parnam{imodel}{ipar}));
    set(gca,'YTick',0.6:0.1:0.9);
    ylabel('fraction correct');
    helper_ILL(gcf,pbar,HGT,false);
    
    % plot inset with coefficient estimation
    figname = sprintf('accuracy_multiple_regression_model_%s_%s_coeff',modnam{imodel},parnam{imodel}{ipar});
    figure('Color','white','Name',figname)
    pbar = 0.5; HGT = 1.5;
    hold on
    switch parnam{imodel}{ipar}
        case 'h'
            icoeff = 3; ylim([-2,0]); yt = [-2,-1,0]; ytl = {'-2' '-1' '0'};
        case 'siginf'
            icoeff = 4; ylim([-0.06,0]); yt = -0.06:0.03:0; ytl = {'-0.06' '-0.03','0'};
        case {'plaps'}
            icoeff = [5 6]; ylim([-0.25,0]); yt = [-.2 0]; ytl = {'-0.2' '0'};
    end
    bavg = nlm.Coefficients.Estimate(icoeff);
    berr = nlm.Coefficients.SE(icoeff);
    if numel(icoeff) < 2
        xlim([0.6,1.4]);
        bar(1,bavg,0.5,'EdgeColor',rgb(imodel,:),'FaceColor',0.5*(rgb(imodel,:)+1),'LineWidth',1);
        plot([1,1],bavg+berr*[-1,+1],'k-');
        xlab = sprintf('%s',nlm.CoefficientNames{icoeff-1});
    else
        pbar = 2*pbar;
        xlim([0.6,2.4]);
        xlab = '';
        for ic = 1:numel(icoeff)
            bar(ic,bavg(ic),0.5,'EdgeColor',rgb(imodel,:),'FaceColor',0.5*(rgb(imodel,:)+1),'LineWidth',1);
            plot([ic,ic],bavg(ic)+berr(ic)*[-1,+1],'k-');
            xlab = sprintf('%s %s',xlab,nlm.CoefficientNames{icoeff(ic)-1});
        end
    end
    set(gca,'XAxisLocation','top','XTick',[]);
    set(gca,'YTick',yt,'YTickLabel',ytl);
    xlabel(sprintf('%s',xlab));
    helper_ILL(gcf,pbar,HGT,false);
end
fprintf('\n')

%% SUPPLEMENTARY: COMPARISON BETWEEN FLAT UNCONDITIONAL, HIERARCHICAL WEIGHTED AND HIERARCHICAL CONDITIONAL INFERENCE MODELS (Supplementary Figure 9)
fprintf('<strong>*** SUPPLEMENTARY: COMPARISON FLAT UNCONDITIONAL, HIERARCHICAL WEIGHTED AND HIERARCHICAL CONDITIONAL INFERENCE MODELS ***</strong>\n\n');

clearvars;
savepath = '../results';

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapse
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapse
    ]/255;
rgb_flatinf = [128,128,128]/255;
rgb_infwght = [150 200 100]/255;
imodel = 3; % infdisc
rgb_infdisc = rgb(imodel,:);

% load fitting output FLAT UNCONDITIONAL INFERENCE model
load(fullfile(savepath,'expe1_results_fit_flat_inference_model.mat'));
out_fit_flatinf1 = out_fit;
ndata1 = size(out_fit_flatinf1,1);
load(fullfile(savepath,'expe2_results_fit_flat_inference_model.mat'));
out_fit_flatinf2 = out_fit;
out_fit_flatinf = cat(1,out_fit_flatinf1,out_fit_flatinf2);
% load fitting output HIERARCHICAL WEIGHTED INFERENCE from stabilizing
% models (EXPE1+EXPE2)
load(fullfile(savepath,'expe1_results_fit_weighted_inference_model.mat'),'out_fit');
out_fit_infwght1 = out_fit;
load(fullfile(savepath,'expe2_results_fit_weighted_inference_model.mat'),'out_fit');
out_fit_infwght2 = out_fit;
out_fit_infwght = cat(1,out_fit_infwght1,out_fit_infwght2);
% load fitting output HIERARCHICAL CONDITIONAL INFERENCE from stabilizing
% models (EXPE1+EXPE2)
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');
out_fit_infdisc1 = out_fit(:,imodel);
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');
out_fit_infdisc2 = out_fit(:,imodel);
out_fit_infdisc = cat(1,out_fit_infdisc1,out_fit_infdisc2);

out_fit = cat(2,out_fit_flatinf,out_fit_infwght,out_fit_infdisc);
ndata = size(out_fit,1);

% BAYESIAN MODEL SELECTION
elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
cfg = [];
cfg.figname = 'expe_merged_bms_flatinf_infwght_infdisc';
cfg.rgb = [rgb_flatinf;rgb_infwght;rgb(imodel,:)];
cfg.elbo = elbo;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
helper_plot_bms(cfg);

% BEST-FITTING MODEL SIMULATIONS
% participants
rep_dat1 = cell2mat(cellfun(@(x)x.crep_sub,out_fit_infdisc(1:ndata1,:),'UniformOutput',0)');
rev_dat1 = cell2mat(cellfun(@(x)x.crev_sub,out_fit_infdisc(1:ndata1,:),'UniformOutput',0)');
rep_dat2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crep_w_sub,out_fit_infdisc((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rev_dat2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crev_w_sub,out_fit_infdisc((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rep_dat = cat(2,rep_dat1,rep_dat2);
rev_dat = cat(2,rev_dat1,rev_dat2);
% flatinf
rep_avg1 = cell2mat(cellfun(@(x)x.crep_avg,out_fit_flatinf(1:ndata1,:),'UniformOutput',0)');
rev_avg1 = cell2mat(cellfun(@(x)x.crev_avg,out_fit_flatinf(1:ndata1,:),'UniformOutput',0)');
rep_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crep_w_avg,out_fit_flatinf((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rev_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crev_w_avg,out_fit_flatinf((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rep_flatinf = cat(2,rep_avg1,rep_avg2);
rev_flatinf = cat(2,rev_avg1,rev_avg2);
% infwght
rep_avg1 = cell2mat(cellfun(@(x)x.crep_avg,out_fit_infwght(1:ndata1,:),'UniformOutput',0)');
rev_avg1 = cell2mat(cellfun(@(x)x.crev_avg,out_fit_infwght(1:ndata1,:),'UniformOutput',0)');
rep_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crep_w_avg,out_fit_infwght((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rev_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crev_w_avg,out_fit_infwght((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rep_infwght = cat(2,rep_avg1,rep_avg2);
rev_infwght = cat(2,rev_avg1,rev_avg2);
% infdisc
rep_avg1 = cell2mat(cellfun(@(x)x.crep_avg,out_fit_infdisc(1:ndata1,:),'UniformOutput',0)');
rev_avg1 = cell2mat(cellfun(@(x)x.crev_avg,out_fit_infdisc(1:ndata1,:),'UniformOutput',0)');
rep_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crep_w_avg,out_fit_infdisc((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rev_avg2 = cell2mat(cellfun(@(x)mean(x,2),cellfun(@(x)x.crev_w_avg,out_fit_infdisc((ndata1+1):end,:),'UniformOutput',0),'UniformOutput',0)');
rep_infdisc = cat(2,rep_avg1,rep_avg2);
rev_infdisc = cat(2,rev_avg1,rev_avg2);

% plot switch curves for the three models
cfg = [];
cfg.figname = 'expe_merged_model_flatinf_SWI_bestfit';
cfg.rgb = rgb_flatinf;
cfg.yavg = 1-mean(rep_flatinf,2);
cfg.yerr = std(rep_flatinf,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_flatinf;
cfg.Yavg = 1-mean(rep_dat,2);
cfg.Yerr = std(rep_dat,[],2)/sqrt(ndata);
helper_plot_rep(cfg)
%
cfg = [];
cfg.figname = 'expe_merged_model_infwght_SWI_bestfit';
cfg.rgb = rgb_infwght;
cfg.yavg = 1-mean(rep_infwght,2);
cfg.yerr = std(rep_infwght,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_infwght;
cfg.Yavg = 1-mean(rep_dat,2);
cfg.Yerr = std(rep_dat,[],2)/sqrt(ndata);
helper_plot_rep(cfg)
%
cfg = [];
cfg.figname = 'expe_merged_model_infdisc_SWI_bestfit';
cfg.rgb = rgb_infdisc;
cfg.yavg = 1-mean(rep_infdisc,2);
cfg.yerr = std(rep_infdisc,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_infdisc;
cfg.Yavg = 1-mean(rep_dat,2);
cfg.Yerr = std(rep_dat,[],2)/sqrt(ndata);
helper_plot_rep(cfg)
% plot reversal curves for the three models
cfg = [];
cfg.figname = 'expe_merged_model_flatinf_REV_bestfit';
cfg.rgb = rgb_flatinf;
cfg.yavg = mean(rev_flatinf,2);
cfg.yerr = std(rev_flatinf,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_flatinf;
cfg.Yavg = mean(rev_dat,2);
cfg.Yerr = std(rev_dat,[],2)/sqrt(ndata);
helper_plot_rev(cfg)
%
cfg = [];
cfg.figname = 'expe_merged_model_infwght_REV_bestfit';
cfg.rgb = rgb_infwght;
cfg.yavg = mean(rev_infwght,2);
cfg.yerr = std(rev_infwght,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_infwght;
cfg.Yavg = mean(rev_dat,2);
cfg.Yerr = std(rev_dat,[],2)/sqrt(ndata);
helper_plot_rev(cfg)
%
cfg = [];
cfg.figname = 'expe_merged_model_infdisc_REV_bestfit';
cfg.rgb = rgb_infdisc;
cfg.yavg = mean(rev_infdisc,2);
cfg.yerr = std(rev_infdisc,[],2)/sqrt(ndata);
cfg.rgb_dat = rgb_infdisc;
cfg.Yavg = mean(rev_dat,2);
cfg.Yerr = std(rev_dat,[],2)/sqrt(ndata);
helper_plot_rev(cfg)

% COMPARE BEST-FITTING PARAMETERS
fprintf('Compare best-fitting parameters (hierarchical weighted versus hierarchical conditional inference\n')
parnam = out_fit_infwght{1}.xnam;
npar = numel(parnam);
xavg_infwght = cell2mat(cellfun(@(s)getfield(s,'xavg'),out_fit_infwght,'UniformOutput',0));
xstd_infwght = cell2mat(cellfun(@(s)getfield(s,'xstd'),out_fit_infwght,'UniformOutput',0));
xavg_infdisc = cell2mat(cellfun(@(s)getfield(s,'xavg'),out_fit_infdisc,'UniformOutput',0));
xstd_infdisc = cell2mat(cellfun(@(s)getfield(s,'xstd'),out_fit_infdisc,'UniformOutput',0));
for ipar = 1:npar
    cfg = [];
    cfg.figname = sprintf('expe_merged_model_infwght_infdisc_CORR_%s',parnam{ipar});
    cfg.xavg = xavg_infwght(:,ipar);
    cfg.xstd = xstd_infwght(:,ipar);
    cfg.yavg = xavg_infdisc(:,ipar);
    cfg.ystd = xstd_infdisc(:,ipar);
    cfg.rgb = mean([rgb_infdisc;rgb_infwght]);
    % xlim/xlabel
    if strcmp(parnam{ipar},'h')
        cfg.xlim = [0 0.4]; cfg.xtick = 0:.1:.4;
        cfg.ylim = [0 0.4]; cfg.ytick = 0:.1:.4;
    elseif strcmp(parnam{ipar},'siginf')
        cfg.xlim = [0 4]; cfg.xtick = 0:1:4; cfg.xticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
        cfg.ylim = [0 4]; cfg.ytick = 0:1:4; cfg.yticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(parnam{ipar},'delta')
        cfg.xlim = [0 3]; cfg.xtick = 0:1:3; cfg.xticklabel = {'0' '1.0' '2.0' '3.0'};
        cfg.ylim = [0 3]; cfg.ytick = 0:1:3; cfg.yticklabel = {'0' '1.0' '2.0' '3.0'};
    end
    cfg.xlabel = sprintf('%s',parnam{ipar});
    cfg.ylabel = sprintf('%s',parnam{ipar});
    cfg.identity = true;
    helper_plot_corr(cfg);
    [~,~,RLO,RUP]=corrcoef(cfg.xavg ,cfg.yavg);
    fprintf('\t\tCI = [%+.3f %+.3f]\n\n',RLO(1,2),RUP(1,2));
end

%% SUPPLEMENTARY: COMPARISON BETWEEN WEIGHTED INFERENCE AND CONDITIONAL INFERENCE MODEL EFFECTS (Supplementary Figure 10)
fprintf('<strong>*** SUPPLEMENTARY: COMPARISON BETWEEN WEIGHTED INFERENCE AND CONDITIONAL INFERENCE MODEL EFFECTS ***</strong>\n\n');
clearvars
savepath = '../results';

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapse
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapse
    ]/255;
rgb_infwght = [150 200 100]/255;
rgb_mod = [rgb_infwght;rgb(3,:)];

if ~exist(fullfile(savepath,'expe_merged_predicted_effects_accuracy_infwght_vs_infdisc.mat'),'file')
    
    % load suboptimal (without stabilization) model fitting results from both experiments
    load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit')
    out_fit_nonstab_expe1 = out_fit;
    ndata1 = size(out_fit_nonstab_expe1,1);
    load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit')
    ndata2 = size(out_fit,1);
    out_fit_nonstab = cat(1,out_fit_nonstab_expe1,out_fit); clearvars out_fit;
    % load weighted inference model fitting results from both experiments
    load(fullfile(savepath,'expe1_results_fit_weighted_inference_model.mat'),'out_fit')
    out_fit_infwght1 = out_fit;
    load(fullfile(savepath,'expe2_results_fit_weighted_inference_model.mat'),'out_fit')
    out_fit_infwght2 = out_fit;
    out_fit_infwght = cat(1,out_fit_infwght1,out_fit_infwght2);
    % load conditional inference model fitting results from both experiments
    load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');
    out_fit_infdisc1 = out_fit(:,3);
    load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');
    out_fit_infdisc2 = out_fit(:,3);
    out_fit_infdisc = cat(1,out_fit_infdisc1,out_fit_infdisc2);
    out_fit = cat(2,out_fit_infwght,out_fit_infdisc);
    ndata = size(out_fit,1);
    
    modnam = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',false),'UniformOutput',false);
    nmodel = numel(modnam);
    parnam = out_fit{1}.xnam;
    npar = numel(parnam);
    
    % compute participants overall accuracy and switch rate
    fprintf('Computing participants overall accuracy\n\t')
    pcor_dat_avg = nan(ndata,1);
    pcor_dat_std = nan(ndata,1);
    prep_dat_avg = nan(ndata,1);
    prep_dat_std = nan(ndata,1);
    for idata = 1:ndata
        fprintf('.')
        dat = out_fit{idata,1}.cfg;
        pcor_dat_avg(idata) = mean(dat.r == dat.s);
        pcor_dat_std(idata) = sqrt(pcor_dat_avg(idata)*(1-pcor_dat_avg(idata))/numel(dat.r));
        rep = [nan;dat.r(2:end) == dat.r(1:end-1)];
        rep(isnan(dat.r),:) = 0;
        rep(dat.t == 1,:) = nan;
        prep_dat_avg(idata) = nanmean(rep,1);
        prep_dat_std(idata) = sqrt(prep_dat_avg(idata)*(1-prep_dat_avg(idata))/numel(rep));
    end
    fprintf('done.\n\n')
    
    fprintf('Computing simulated effects of belief stabilization strategies on decision accuracy (long)\n')
    % simulation configuration
    npar_sim = 21;  % number of values taken by stabilization parameter
    nsmp = 1e3;     % number of sample per simulation
    
    % simulation output
    pcor_sim_avg = nan(ndata,npar_sim,nmodel);      % simulated mean accuracy varying delta
    prep_sim_avg = nan(ndata,npar_sim,nmodel);      % simulated mean accuracy varying delta
    pswi_sim_avg = nan(ndata,npar_sim,nmodel,2);    % simulated fraction switch varying delta for correct and error trials
    lastbef_sim_avg = nan(ndata,nmodel);            % simulated mean accuracy on last trial before a reversal varying delta
    firstaf_sim_avg = nan(ndata,nmodel);            % simulated mean accuracy on first trial after a reversal varying delta
    
    % take mean of non-stabilized best-fitting hazard rate and inference noise
    h0      = mean(cellfun(@(s)getfield(s,'h'),out_fit_nonstab));
    siginf0 = mean(cellfun(@(s)getfield(s,'siginf'),out_fit_nonstab));
    
    delta_vec = linspace(0,3,npar_sim);  % parameter values for simulation
    
    for imodel = 1:nmodel
        fprintf('\t%s \t..',modnam{imodel});
        % run simulations
        for idata = 1:ndata
            fprintf('.');
            % get model predictions varying delta
            cfg = out_fit{idata}.cfg;
            cfg_sim = [];
            cfg_sim.t = cfg.t;
            cfg_sim.s = cfg.s;
            cfg_sim.pfal = cfg.pfal;
            if idata > ndata1 % second experiment
                cfg_sim.w = cfg.w;
            end
            cfg_sim.modtype = modnam{imodel};
            cfg_sim.h = h0;
            cfg_sim.siginf = siginf0;
            cfg_sim.sigsel = 0;
            cfg_sim.nsmp = nsmp;
            cfg_sim.seed = idata*100;
            
            cfg_tmp = cfg_sim;
            for isim = 1:npar_sim
                cfg_tmp.delta = delta_vec(isim);
                if idata <= ndata1
                    out = sim_model(cfg_tmp);
                else
                    out = sim_model_expe2(cfg_tmp);
                end
                % get fraction correct
                pcor_sim_avg(idata,isim,imodel) = out.pcor;
                % get fraction switch
                prep_sim_avg(idata,isim,imodel) = out.prep;
                % identify switch trials
                swi_sim = [false(1,nsmp);out.r(2:end,:) ~= out.r(1:end-1,:)];
                swi_sim(cfg_sim.t == 1,:) = false;
                % identify trials where switching is the correct answer
                swi_cor = [false(1,nsmp);(out.r(1:end-1,:) ~= cfg_tmp.s(2:end))];
                swi_cor(cfg_sim.t == 1,:) = false;
                % identify trials where switching is NOT the correct answer
                swi_inc = [false(1,nsmp);(out.r(1:end-1,:) == cfg_tmp.s(2:end))];
                swi_inc(cfg_sim.t == 1,:) = false;
                % compute fraction switch when switching is correct
                pswi_sim_avg(idata,isim,imodel,1) = mean(sum(swi_sim & swi_cor)./nansum(swi_cor));
                % compute fraction switch when switching is incorrect
                pswi_sim_avg(idata,isim,imodel,2) = mean(sum(swi_sim & swi_inc)./nansum(swi_inc));
                
                lastbef_sim_avg(idata,isim,imodel) = out.c.rev_avg(4);
                firstaf_sim_avg(idata,isim,imodel) = out.c.rev_avg(5);
            end
        end
    fprintf('done.\n')
    end

    save(fullfile(savepath,'expe_merged_predicted_effects_accuracy_infwght_vs_infdisc.mat'),...
        'ndata','modnam','nmodel',...
        'parnam','npar','delta_vec',...
        'pcor_dat_avg','pcor_dat_std','prep_dat_avg','prep_dat_std',...
        'pcor_sim_avg','prep_sim_avg','pswi_sim_avg',...
        'lastbef_sim_avg','firstaf_sim_avg')
else
    load(fullfile(savepath,'expe_merged_predicted_effects_accuracy_infwght_vs_infdisc.mat'))

    pbar = .9; HGT = 4;
    
    % plot predicted accuracy varying delta
    fprintf('Optimal performance achieved: \n')
    for imodel = 1:nmodel
        % interpolate to find maximum
        xdat = delta_vec';
        ydat = mean(pcor_sim_avg(:,:,imodel),1)';
        par_opt     = fminbnd(@(x)-interp1(xdat,ydat,x,'spline'),0,max(xdat));
        pcor_opt    = interp1(xdat,ydat,par_opt,'spline');
        pcor0       = interp1(xdat,ydat,0,'spline');
        
        figname = sprintf('accuracy_simulation_model_%s',modnam{imodel});
        figure('Color','white','Name',figname)
        hold on
        xlim([min(xdat) max(xdat)]);
        ylim([0.5 0.9]);
        plot(par_opt*[1,1],ylim,'k-');
        plot(xdat,ydat,'-','LineWidth',1.5,'Color',rgb_mod(imodel,:));
        plot(par_opt,pcor_opt,'ko','MarkerSize',10,'MarkerFaceColor',rgb_mod(imodel,:));
        
        gain = (pcor_opt-pcor0)*100;
        fprintf('\t%s:\t%4.2f%% (= p0%+4.3f) for delta \t= %4.3f\n',modnam{imodel},pcor_opt*100,gain,par_opt)

        set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        plot(xlim,pcor0*[1 1],'k:')
        xlabel('delta');
        set(gca,'YTick',0.5:0.1:0.9);
        ylabel('fraction correct');
        helper_ILL(gcf,pbar,HGT);
    end
    fprintf('\n')
    
    % plot predicted fraction switch varying delta
    for imodel = 1:nmodel
        xdat = delta_vec';
        ydat = 1-mean(prep_sim_avg(:,:,imodel),1)';
        prep0       = interp1(xdat,ydat,0,'spline');
        
        figname = sprintf('switch_rate_simulation_model_%s',modnam{imodel});
        figure('Color','white','Name',figname)
        hold on
        xlim([min(xdat) max(xdat)]);
        ylim([0 0.3]);
        plot(xdat,ydat,'-','LineWidth',1.5,'Color',rgb_mod(imodel,:));

        set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        plot(xlim,prep0*[1 1],'k:')
        xlabel('delta');
        set(gca,'YTick',0:0.1:0.3);
        ylabel('fraction switch');
        helper_ILL(gcf,pbar,HGT);
    end
    
    % plot predicted accuracy on the last trial before and the first trial
    % after a reversal varying delta
    pbar = .9; HGT = 4;
    for imodel = 1:nmodel
        xvec = delta_vec;
        Y1 = squeeze(mean(firstaf_sim_avg(:,:,imodel)));
        Y2 = 1-squeeze(mean(lastbef_sim_avg(:,:,imodel)));
        pcor0_1 = interp1(xvec,Y1,0,'spline');
        pcor0_2 = interp1(xvec,Y2,0,'spline');
        % last before / first after
        figname = sprintf('accuracy_before_after_model_%s_delta',modnam{imodel});
        figure('Color','white','Name',figname)
        hold on
        xlim([min(xvec) max(xvec)]);
        ylim([0 1]);
        plot(xvec,Y1,'Color',rgb_mod(imodel,:),'LineWidth',1.5);
        plot(xvec,Y2,'Color',rgb_mod(imodel,:),'LineWidth',1.5,'LineStyle',':');
        plot(xlim,pcor0_1*[1 1],'Color',[.8 .8 .8])
        plot(xlim,pcor0_2*[1 1],'Color',[.8 .8 .8])
        set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        xlabel('delta');ylabel('fraction correct')
        helper_ILL(gcf,pbar,HGT,false);
    end
    
    % plot corresponding fraction switch on trials where switch is correct or
    % incorrect
    for imodel = 1:nmodel
        par_opt = fminbnd(@(x)-interp1(delta_vec',mean(pcor_sim_avg(:,:,imodel),1)',x,'spline'),0,max(delta_vec));
        
        xdat = delta_vec;
        ydat = squeeze(mean(pswi_sim_avg(:,:,imodel,:),1));
        pswi_opt = interp1(xdat,ydat,par_opt,'spline');
        pswi0 = interp1(xdat,ydat,0,'spline');
        figname = sprintf('fraction_switch_correct_error_simulation_model_%s',modnam{imodel});
        figure('Color','white','Name',figname)
        pbar = .9; HGT = 4;
        hold on
        xlim([min(xdat) max(xdat)]);
        ylim([0 .6]);
        plot(par_opt*[1,1],ylim,'k-');
        plot(par_opt,pswi_opt(1),'ko','MarkerSize',10,'MarkerFaceColor',rgb_mod(imodel,:));
        plot(par_opt,pswi_opt(2),'ko','MarkerSize',10,'MarkerFaceColor',rgb_mod(imodel,:));
        plot(xdat,ydat(:,1),'Color',rgb_mod(imodel,:),'LineWidth',1.5);
        plot(xdat,ydat(:,2),'Color',rgb_mod(imodel,:),'LineWidth',1.5,'LineStyle',':');
        plot(xlim,pswi0(1)*[1 1],'k:');plot(xlim,pswi0(2)*[1 1],'k:');
        set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        xlabel('delta'); ylabel('fraction switch')
        helper_ILL(gcf,pbar,HGT);
    end
end

%% SUPPLEMENTARY: FIT TO EACH CHOICE (Supplementary Figure 11)
fprintf('<strong>*** SUPPLEMENTARY: FIT TO EACH CHOICE ***</strong>\n\n');

clearvars;
savepath = '../results';

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapse
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapse
    ]/255;

% load output from fit to reversal and switch curves
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');
models = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);
ndata = size(out_fit,1);
nmodel = size(out_fit,2);
elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_stab = cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0);
crep_stab = cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0);
% participants
crev_dat = cell2mat(cellfun(@(s)getfield(s,'crev_sub'),out_fit(:,1),'UniformOutput',0)')';
crep_dat = cell2mat(cellfun(@(s)getfield(s,'crep_sub'),out_fit(:,1),'UniformOutput',0)')';

% load output from fit to each choice
load(fullfile(savepath,'expe1_results_fit_stabilizing_models_each_choice.mat'),'out_fit');
elbo_choice = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_choice = cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0);
crep_choice = cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0);
xavg_choice = cellfun(@(x)x(end),cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',0));
xnam_choice = cellfun(@(x)x(end),cellfun(@(x)getfield(x,'xnam'),out_fit(1,:),'UniformOutput',0));

% BAYESIAN MODEL SELECTION
cfg = [];
cfg.figname = 'expe1_bms_choice';
cfg.rgb = rgb;
cfg.elbo = elbo_choice;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
helper_plot_bms(cfg);

% REVERSAL BEHAVIOR
for imodel = 1:nmodel
    % plot stabilization models versus non stabilized model
    cfg = [];
    cfg.figname = sprintf('expe1_REV_model_%s_choice',models{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.rgb_dashed = cfg.rgb;
    cfg.yavg = mean(cell2mat(crev_choice(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crev_choice(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = mean(crev_dat)';
    cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
    cfg.yavg_dashed = mean(cell2mat(crev_stab(:,imodel)'),2);
    helper_plot_rev(cfg)
    
    cfg = [];
    cfg.figname = sprintf('expe1_SWI_model_%s_choice',models{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.rgb_dashed = cfg.rgb;
    cfg.yavg = 1-mean(cell2mat(crep_choice(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crep_choice(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = 1-mean(crep_dat)';
    cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
    cfg.yavg_dashed = 1-mean(cell2mat(crep_stab(:,imodel)'),2);
    helper_plot_rep(cfg)
    
    fprintf('\t%s',xnam_choice{imodel});
    fprintf('\t%4.3f +/- %4.3f\n',mean(xavg_choice(:,imodel)),std(xavg_choice(:,imodel),[],1)/sqrt(ndata));
end
