%% EXPERIMENT 1 RUNNER ANALYSIS STABILIZATION
% needed on the path for running analysis:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)

% clear workspace
clearvars;
addpath '../helper_plot';
addpath '../matlab_functions';

savepath = '../results/';

% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapses
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapses
    ]/255;

%% PREDICTED EFFECTS OF RESPONSE STABILIZATION STRATEGIES (Figure 3)
if exist(fullfile(savepath,'expe1_predicted_effects_stabilization.mat'),'file')
    load(fullfile(savepath,'expe1_predicted_effects_stabilization.mat'))
else
    load(fullfile(savepath,'expe1_results_fit_stabilizing_parameters.mat'))
    ndata = size(out_fit,1);
    modnam = cellfun(@(s)getfield(s,'modtype'),out_fit(1,:),'UniformOutput',0);
    nmodel = numel(modnam);
    parnam = cellfun(@(c)c{3},cellfun(@(s)getfield(s,'xnam'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);
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
    load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit');
    % response switch curves and response reversal curves
    crev_stab_ex_ante = cell(ndata,nmodel);
    crep_stab_ex_ante = cell(ndata,nmodel);
    for imodel = 1:nmodel
        fprintf('\t%s \t..',models{imodel});
        for idata = 1:ndata
            fprintf('.')
            cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
            cfg_sim = [];
            cfg_sim.t = cfg.t;
            cfg_sim.s = cfg.s;
            cfg_sim.h = out_fit{idata}.h;
            cfg_sim.siginf = out_fit{idata}.siginf;
            cfg_sim.modtype = modnam{imodel};
            cfg_sim.sigsel = 0;
            cfg_sim.pfal = 0.20;
            cfg_sim.nsmp = 1e3;
            % add stabilizing parameter
            cfg_sim.(parnam{imodel}) = stab_params(idata,imodel);
            out_sim = sim_model(cfg_sim);
            crev_stab_ex_ante{idata,imodel} = out_sim.c.rev_avg;
            crep_stab_ex_ante{idata,imodel} = out_sim.c.rep_avg;
        end
        fprintf('done.\n')
    end
    fprintf('\n');
    % Compute simulated effects of the 5 candidate strategies on response 
    % switches just before and after a reversal
    npar_sim = 20; % number of sample of stabilizing parameter to simulate
    stab_param_max = [0.5 .7 2 1.5 .7];   % max stab parameter values for simulation
    stab_param_vec = nan(nmodel,npar_sim); % parameter values for simulation
    last_before_avg = nan(nmodel,npar_sim);
    first_after_avg = nan(nmodel,npar_sim);
    last_before_err = nan(nmodel,npar_sim);
    first_after_err = nan(nmodel,npar_sim);
    last_before_ex_ante = nan(ndata,nmodel);
    first_after_ex_ante = nan(ndata,nmodel);
    fprintf('<strong>Computing fraction switch on last trial before and first trial after each reversal (long)</strong>\n')
    for imodel = 1:nmodel
        stab_param_vec(imodel,:) = linspace(0,stab_param_max(imodel),npar_sim);
        fprintf('\t%s \t..',models{imodel});
        for isim = 1:npar_sim
            fprintf('.')
            last_before_tmp = nan(1,ndata);
            first_after_tmp = nan(1,ndata);
            for idata = 1:ndata
                cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
                cfg_sim = [];
                cfg_sim.t = cfg.t;
                cfg_sim.s = cfg.s;
                cfg_sim.h = out_fit{idata}.h;
                cfg_sim.siginf = out_fit{idata}.siginf;
                cfg_sim.modtype = modnam{imodel};
                cfg_sim.sigsel = 0;
                cfg_sim.pfal = 0.20;
                cfg_sim.nsmp = 1e3;
                % add stabilizing parameter
                cfg_sim.(parnam{imodel}) = stab_param_vec(imodel,isim);
                out_sim = sim_model(cfg_sim);
                % extract fraction switch
                last_before_tmp(idata) = 1-out_sim.c.rep_avg(4);
                first_after_tmp(idata) = 1-out_sim.c.rep_avg(5);
            end
            last_before_avg(imodel,isim) = mean(last_before_tmp);
            first_after_avg(imodel,isim) = mean(first_after_tmp);
            last_before_err(imodel,isim) = std(last_before_tmp)/sqrt(ndata);
            first_after_err(imodel,isim) = std(first_after_tmp)/sqrt(ndata);
        end
        fprintf('done\n')
        % predictions for mean stabilization parameters matching overall
        % accuracy
        for idata = 1:ndata
            cfg = out_fit{idata}.cfg; % suboptimal model (without stabilization)
            cfg_sim = [];
            cfg_sim.t = cfg.t;
            cfg_sim.s = cfg.s;
            cfg_sim.h = out_fit{idata}.h;
            cfg_sim.siginf = out_fit{idata}.siginf;
            cfg_sim.modtype = modnam{imodel};
            cfg_sim.sigsel = 0;
            cfg_sim.pfal = 0.20;
            cfg_sim.nsmp = 1e3;
            % add stabilizing parameter
            cfg_sim.(parnam{imodel}) = mean(stab_params(:,imodel));
            out_sim = sim_model(cfg_sim);
            % extract fraction switch
            last_before_ex_ante(idata,imodel) = 1-out_sim.c.rep_avg(4);
            first_after_ex_ante(idata,imodel) = 1-out_sim.c.rep_avg(5);
        end
    end
    fprintf('\n');
    save(fullfile(savepath,'expe1_predicted_effects_stabilization.mat'),...
        'ndata','modnam','parnam','nmodel','models','stab_params',...
        'crev_stab_ex_ante','crep_stab_ex_ante','last_before_ex_ante','first_after_ex_ante',...
        'stab_param_vec',...
        'last_before_avg','first_after_avg','last_before_err','first_after_err')
end

% extract model fitting output (without stabilization) to be stabilized
load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit');
crev_nonstab = cell2mat(cellfun(@(c)c',cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0),'UniformOutput',0));
crep_nonstab = cell2mat(cellfun(@(c)c',cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0),'UniformOutput',0));

% Plot simulated effects of the 5 candidate strategies on reversal behavior
% (Figure 3a) 
fprintf('<strong>Stabilization parameter values best-fitting participants'' overall switch rate:</strong>\n');

for imodel = 1:nmodel
    % plot stabilization models versus non stabilized model
    cfg = [];
    cfg.figname = sprintf('Figure3a_expe1_ex_ante_model_%s_SWI',modnam{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.yavg = 1-mean(cell2mat(crep_stab_ex_ante(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crep_stab_ex_ante(:,imodel)'),[],2)/sqrt(ndata);
    cfg.yavg_dashed = mean(1-crep_nonstab,1);
    helper_plot_rep(cfg)
    
    cfg = [];
    cfg.figname = sprintf('Figure3a_expe1_ex_ante_model_%s_REV',modnam{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.yavg = mean(cell2mat(crev_stab_ex_ante(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crev_stab_ex_ante(:,imodel)'),[],2)/sqrt(ndata);
    cfg.yavg_dashed = mean(crev_nonstab,1);
    helper_plot_rev(cfg)

    fprintf('\t%s \t%s\t= %5.3f +/- %5.3f\n' ,models{imodel},parnam{imodel},mean(stab_params(:,imodel)),std(stab_params(:,imodel))/sqrt(ndata));
end
fprintf('\n')

% Plot simulated effects of the 5 candidate strategies on response switches
% just before and after a reversal (Figure 3b)
pbar = 1.1;
HGT = 4;
for imodel = 1:nmodel
    % plot reversal curve
    figname = sprintf('Figure3b_expe1_ex_ante_model_%s_firstlast',modnam{imodel});
    figure('Color','white','Name',figname);
    
    Y1 = last_before_avg(imodel,:);
    Y2 = first_after_avg(imodel,:);
    Y1err = last_before_err(imodel,:);
    Y2err = first_after_err(imodel,:);
    xvec = stab_param_vec(imodel,:);
    xl = [min(xvec) max(xvec)];
    
    Y1_ex_ante = mean(last_before_ex_ante(:,imodel));
    Y2_ex_ante = mean(first_after_ex_ante(:,imodel));
    
    hold on
    xlim(xl);
    ylim([0,.4]);
    patch([xvec,fliplr(xvec)],[Y1+Y1err,fliplr(Y1-Y1err)], ...
        0.5*(rgb(imodel,:)+1),'EdgeColor','none');
    plot(xvec,Y1,'Color',rgb(imodel,:),'LineWidth',1.5,'LineStyle',':');
    patch([xvec,fliplr(xvec)],[Y2+Y2err,fliplr(Y2-Y2err)], ...
        0.5*(rgb(imodel,:)+1),'EdgeColor','none');
    plot(xvec,Y2,'Color',rgb(imodel,:),'LineWidth',1.5);
    plot([mean(stab_params(:,imodel)),mean(stab_params(:,imodel))],ylim,'k-'); % mean parameter value
    plot(mean(stab_params(:,imodel)),Y1_ex_ante,'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
    plot(mean(stab_params(:,imodel)),Y2_ex_ante,'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
    set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
    set(gca,'YTick',0:0.1:.4);
    xlabel(sprintf('%s',parnam{imodel}));
    ylabel('fraction switch');
    helper_ILL(gcf,pbar,HGT);
end    

%% CONFUSION MATRIX for EX-ANTE MODEL RECOVERY (Figure 4a)
load(fullfile(savepath,'expe1_results_recovery_stabilizing_models_elbo.mat'))
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
figure('Color','white','Name','Figure4a_expe1_confusion_matrix_stabilizing_models');
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

%% STABILIZING MODELS - extract fitting output
load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit');
modnam = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);
ndata = size(out_fit,1);
nmodel = size(out_fit,2);

elbo = cellfun(@(s)getfield(s,'elbo'),out_fit);
crev_stab = cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0);
crep_stab = cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0);

% stabilization parameters value
xavg = cellfun(@(x)x(3),cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',0));
xnam = cellfun(@(x)x(3),cellfun(@(x)getfield(x,'xnam'),out_fit(1,:),'UniformOutput',0));

% extract model fitting output (without stabilization)
load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit');
crev_nonstab = cell2mat(cellfun(@(s)getfield(s,'crev_avg'),out_fit,'UniformOutput',0)')';
crep_nonstab = cell2mat(cellfun(@(s)getfield(s,'crep_avg'),out_fit,'UniformOutput',0)')';

% participants
crev_dat = cell2mat(cellfun(@(s)getfield(s,'crev_sub'),out_fit,'UniformOutput',0)')';
crep_dat = cell2mat(cellfun(@(s)getfield(s,'crep_sub'),out_fit,'UniformOutput',0)')';

%% BAYESIAN MODEL SELECTION (Figure 4b)
cfg = [];
cfg.figname = 'Figure_4b_expe1_bms';
cfg.rgb = rgb;
cfg.elbo = elbo;
cfg.ylabel = 'model probability';
cfg.xlabel = 'model';
helper_plot_bms(cfg);

%% REVERSAL BEHAVIOR - BEST-FITTING STABILIZING MODEL SIMULATIONS (Figure 4c)
fprintf('<strong>Stabilization parameter values best-fitting participants'' reversal behavioral metrics:</strong>\n')
for imodel = 1:nmodel
    % plot stabilization models versus non stabilized model    
    cfg = [];
    cfg.figname = sprintf('Figure_4c_expe1_SWI_model_%s',modnam{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = 1-mean(cell2mat(crep_stab(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crep_stab(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = 1-mean(crep_dat)';
    cfg.Yerr = std(crep_dat,[],1)'/sqrt(ndata);
    cfg.yavg_dashed = 1-mean(crep_nonstab)';
    helper_plot_rep(cfg)
    
    cfg = [];
    cfg.figname = sprintf('Figure_4c_expe1_REV_model_%s',modnam{imodel});
    cfg.rgb = rgb(imodel,:);
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = mean(cell2mat(crev_stab(:,imodel)'),2);
    cfg.yerr = std(cell2mat(crev_stab(:,imodel)'),[],2)/sqrt(ndata);
    cfg.Yavg = mean(crev_dat)';
    cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
    cfg.yavg_dashed = mean(crev_nonstab)';
    helper_plot_rev(cfg)

    fprintf('\t%s',xnam{imodel});
    fprintf('\t%4.3f +/- %4.3f\n',mean(xavg(:,imodel)),std(xavg(:,imodel),[],1)/sqrt(ndata));
end
fprintf('\n')
