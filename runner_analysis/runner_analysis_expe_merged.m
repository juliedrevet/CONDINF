%% RUNNER ANALYSIS EXPERIMENTS MERGED
addpath '../helper_plot';
addpath '../matlab_functions';

fprintf('<strong>***\tEXPERIMENTS 1 and 2\t***</strong>\n\n')
%% INTERINDIVIDUAL VARIABILITY IN CONDITIONAL INFERENCE (Figure 7)
clearvars;
savepath = '../results';

fprintf('<strong>Interindividual variability in conditional inference</strong>\n\n')
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

pbar = 1.1; HGT = 4;
fprintf('Best-fitting model parameter values:\n')
for ipar = 1:npar
    X = xavg_fit(:,ipar);
    
    figname = sprintf('Figure_7a_conditional_inference_parameter_%s',parnam{ipar});
    figure('Color','white','Name',figname);
    % plot histogram of best-fitting parameters
    h_par = histogram(X,'NumBins',8,'EdgeColor','w','FaceColor',(1+rgb(imodel,:))/2,'Normalization','pdf');
    hold on
    xk = linspace(min(0),max(X),100);
    pk = ksdensity(X,xk);
    plot(xk,pk,'Color',rgb(imodel,:),'LineWidth',2)
    
    line(median(X)*[1 1],[0 max(ylim)],'LineStyle',':','Color','k')
    text(median(X),.9*max(ylim),sprintf('%4.3f',median(X)));
    xlabel(sprintf('%s',parnam{ipar}))
    
    if strcmp(parnam{ipar},'h')
        set(gca,'XTick',0:.1:.4,'XTickLabel',{'0' '0.1' '0.2' '0.3' '0.4'});
    elseif strcmp(parnam{ipar},'siginf')
        set(gca,'XTick',0:1.0:4.0,'XTickLabel',{'0' '1.0' '2.0' '3.0' '4.0'});
    elseif strcmp(parnam{ipar},'delta')
        set(gca,'XTick',0:1.0:3.0,'XTickLabel',{'0' '1.0' '2.0' '3.0'});
    end
    
    fprintf('\tmedian %s\t= %4.5f\n',parnam{ipar},median(X))
    
    helper_ILL(gcf,pbar,HGT,4);
end

fprintf('\nBest-fitting model parameter correlations:\n')
for ipar = 1:npar
    X = xavg_fit(:,ipar);
    % plot correlation between conditional inference model parameters
    cfg = [];
    cfg.figname = sprintf('Figure_7b_conditional_inference_corr_%s_%s',parnam{ipar},parnam{mod(ipar,3)+1});
    cfg.rgb = rgb(imodel,:);
    cfg.xavg = X;
    cfg.yavg = xavg_fit(:,mod(ipar,3)+1);
    cfg.ystd = xstd_fit(:,mod(ipar,3)+1);
    % xlim/xlabel
    if strcmp(parnam{ipar},'h')
        cfg.xlim = [0 0.4]; cfg.xtick = 0:.1:.4;
    elseif strcmp(parnam{ipar},'siginf')
        cfg.xlim = [0 4]; cfg.xtick = 0:1:4; cfg.xticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(parnam{ipar},'delta')
        cfg.xlim = [0 3]; cfg.xtick = 0:1:3; cfg.xticklabel = {'0' '1.0' '2.0' '3.0'};
    end
    % ylim/ylabel
    if strcmp(parnam{mod(ipar,3)+1},'h')
        cfg.ylim = [0 .4]; cfg.ytick = 0:.1:.4;
    elseif strcmp(parnam{mod(ipar,3)+1},'siginf')
        cfg.ylim = [0 4]; cfg.ytick = 0:1:4; cfg.yticklabel = {'0' '1.0' '2.0' '3.0' '4.0'};
    elseif strcmp(parnam{mod(ipar,3)+1},'delta')
        cfg.ylim = [0 3]; cfg.ytick = 0:1:3; cfg.yticklabel = {'0' '1.0' '2.0' '3.0'};
    end
    cfg.xlabel = sprintf('%s',parnam{ipar});
    cfg.ylabel = sprintf('%s',parnam{mod(ipar,3)+1});
    cfg.identity = false;
    helper_plot_corr(cfg);
    
    [R,P,RLO,RUP]=corrcoef(cfg.xavg ,cfg.yavg);
    fprintf('\t\tCI = [%+4.5f %+4.5f]\n',RLO(1,2),RUP(1,2))
end


for ipar = 1:npar
    % mediansplit on best-fitting parameter values
    X = xavg_fit(:,ipar);
    xmed = median(X);
    idx = X <= xmed;
    
    % RESPONSE SWITCH curves
    crep_mod_avg = cat(2,mean(crep_mod(idx,:))',mean(crep_mod(~idx,:))');
    crep_mod_err = cat(2,serr(crep_mod(idx,:),1)',serr(crep_mod(~idx,:),1)');
    crep_dat_avg = cat(2,mean(crep_dat(idx,:))',mean(crep_dat(~idx,:))');
    crep_dat_err = cat(2,serr(crep_dat(idx,:),1)',serr(crep_dat(~idx,:),1)');
    cfg = [];
    cfg.figname = sprintf('Figure_7c_conditional_inference_mediansplit_SWI_%s',parnam{ipar});
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
    cfg.figname = sprintf('Figure_7c_conditional_inference_mediansplit_REV_%s',parnam{ipar});
    cfg.rgb = [(1+rgb(imodel,:))/2;rgb(imodel,:)];
    cfg.rgb_dat = cfg.rgb;
    cfg.yavg = crev_mod_avg;
    cfg.yerr = crev_mod_err;
    cfg.Yavg = crev_dat_avg;
    cfg.Yerr = crev_dat_err;
    helper_plot_rev(cfg)
end
% variance explained
fprintf('\n<strong>Variance explained by conditional inference model</strong>\n');
xdat = cat(2,crev_dat,crep_dat); % human data
vexpl_dat = nan(3,2); % variance explained for reversal=1 and repetition=2 curves for subject data
xres = nan(size(xdat));
modelfun = @(b,x)b(1)+b(2)*x; % model function
start = [0;0]; % starting values
for ipar = 1:3 % lignes = de vexpl = quel paramètre
    X = xavg_fit(:,ipar);
    for j = 1:16
        Y = xdat(:,j);
        nlm = fitnlm(X,Y,modelfun,start);
        xres(:,j) = nlm.Residuals.Raw;
    end
    vexpl_dat(ipar,1) = 1-sum(var(xres(:,1:8)))/sum(var(xdat(:,1:8)));
    vexpl_dat(ipar,2) = 1-sum(var(xres(:,9:16)))/sum(var(xdat(:,9:16)));
    fprintf('Parameter %s:\n',parnam{ipar})
    fprintf('\t* switch curve:   %3.1f%%\n',vexpl_dat(ipar,2)*100);
    fprintf('\t* reversal curve: %3.1f%%\n\n',vexpl_dat(ipar,1)*100);
end

%% INCREASED DECISION ACCURACY THROUGH CONDITIONAL INFERENCE (Figure 8)
clearvars
savepath = '../results';
% define model colors for plots
rgb = [ ...
    238,197,141; ... % senbias: sensory bias
    048,102,158; ... % inflaps: inference lapses
    134,198,226; ... % infdisc: conditional inference
    171,131,171; ... % selbias: repetition bias
    099,074,135  ... % selepsi: response lapses
    ]/255;

fprintf('<strong>Increased decision accuracy through conditional inference</strong>\n\n')
if exist(fullfile(savepath,'expe_merged_predicted_effects_accuracy.mat'),'file')
    load(fullfile(savepath,'expe_merged_predicted_effects_accuracy.mat'))
else
    % load suboptimal (without stabilization) model fitting results from both experiments
    load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'),'out_fit')
    out_fit_nonstab_expe1 = out_fit;
    ndata1 = size(out_fit_nonstab_expe1,1);
    load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit')
    out_fit_nonstab = cat(1,out_fit_nonstab_expe1,out_fit); clearvars out_fit;
    % load stabilizing models fitting results from both experiments
    load(fullfile(savepath,'expe1_results_fit_stabilizing_models.mat'),'out_fit')
    out_fit_expe1 = out_fit;
    load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit')
    ndata2 = size(out_fit,1);
    out_fit = cat(1,out_fit_expe1,out_fit);
    ndata = size(out_fit,1);
    modnam = cellfun(@(s)getfield(s,'modtype'),cellfun(@(s)getfield(s,'cfg'),out_fit(1,:),'UniformOutput',false),'UniformOutput',false);
    nmodel = numel(modnam);
    models = {'Perceptual bias','Inference lapses','Conditional inference','Repetition bias','Response lapses'};
    xavg_fit = cellfun(@(s)getfield(s,'xavg'),out_fit,'UniformOutput',0);
    xstd_fit = cellfun(@(s)getfield(s,'xstd'),out_fit,'UniformOutput',0);
    parnam = cellfun(@(s)getfield(s,'xnam'),out_fit(1,:),'UniformOutput',0);
    stab_param_nam = cellfun(@(c)c{3},cellfun(@(s)getfield(s,'xnam'),out_fit(1,:),'UniformOutput',0),'UniformOutput',0);
    npar = numel(parnam{1});
    % compute participants overall accuracy
    fprintf('Computing participants overall accuracy\n\t')
    pcor_dat_avg = nan(ndata,1);
    pcor_dat_std = nan(ndata,1);
    for idata = 1:ndata
        fprintf('.')
        dat = out_fit{idata,1}.cfg;
        pcor_dat_avg(idata) = mean(dat.r == dat.c);
        pcor_dat_std(idata) = sqrt(pcor_dat_avg(idata)*(1-pcor_dat_avg(idata))/numel(dat.r));
    end
    fprintf('done.\n\n')
    
    fprintf('Computing simulated effects of belief stabilization mechanisms on decision accuracy (long)\n')
    % simulation configuration
    npar_sim = 21;  % number of values taken by stabilization parameter
    nsmp = 1e3;     % number of sample per simulation
    
    % simulation output
    pcor_mod_avg = nan(ndata,nmodel); % predicted best-fitting model mean accuracy
    pcor_mod_std = nan(ndata,nmodel);
    pcor_sim_avg = nan(ndata,npar_sim,nmodel);      % simulated mean accuracy varying stab parameter
    pswi_sim_avg = nan(ndata,npar_sim,nmodel,2);    % simulated fraction switch varying stab parameter for correct and error trials
    lastbef_sim_avg = nan(ndata,nmodel);            % simulated mean accuracy on last trial before a reversal varying stab parameter
    firstaf_sim_avg = nan(ndata,nmodel);            % simulated mean accuracy on first trial after a reversal varying stab parameter
    stab_param_max = [1.5 1 3 3 1];         % max stab parameter values for simulation
    stab_param_vec = nan(nmodel,npar_sim);  % parameter values for simulation
    
    % take mean of non-stabilized best-fitting hazard rate and inference noise
    h0      = mean(cellfun(@(s)getfield(s,'h'),out_fit_nonstab));
    siginf0 = mean(cellfun(@(s)getfield(s,'siginf'),out_fit_nonstab));
    
    for imodel = 1:nmodel
        fprintf('\t%s \t..',models{imodel});
        
        stab_param_vec(imodel,:) = linspace(0,stab_param_max(imodel),npar_sim);
        % run simulations
        for idata = 1:ndata
            fprintf('.');
            % get model predictions from participants fitted values
            cfg_sim      = out_fit{idata,imodel}.cfg;
            cfg_sim.nsmp = nsmp;
            cfg_sim.seed = idata*100;
            for ipar = 1:npar
                cfg_sim.(parnam{imodel}{ipar}) = out_fit{idata,imodel}.(parnam{imodel}{ipar});
            end
            if idata <= ndata1
                out = sim_model(cfg_sim);
            else
                out = sim_model_expe2(cfg_sim);
            end
            pcor_mod_avg(idata,imodel) = out.pcor;
            pcor_mod_std(idata,imodel) = std(mean(out.r == cfg_sim.s));
            % get model predictions varying stabilization parameter
            cfg_sim         = out_fit{idata,imodel}.cfg;
            cfg_sim.nsmp    = nsmp;
            cfg_sim.seed    = idata*100;
            cfg_sim.h       = h0;
            cfg_sim.siginf  = siginf0;
            cfg_tmp = cfg_sim;
            for isim = 1:npar_sim
                cfg_tmp.(stab_param_nam{imodel}) = stab_param_vec(imodel,isim);
                if idata <= ndata1
                    out = sim_model(cfg_tmp);
                else
                    out = sim_model_expe2(cfg_tmp);
                end
                pcor_sim_avg(idata,isim,imodel) = out.pcor;
                % identify switch trials
                swi_sim = [false(1,nsmp);out.r(2:end,:) ~= out.r(1:end-1,:)];
                swi_sim(cfg_sim.t == 1,:) = false;
                % identify trials where switching is the correct answer
                swi_cor = [false(1,nsmp);(out.r(1:end-1,:) ~= cfg_tmp.s(2:end))];
                swi_cor(cfg_sim.t == 1,:) = false;
                % identify trials where switching is NOT the correct answer
                swi_inc = [false(1,nsmp);(out.r(1:end-1,:) == cfg_tmp.c(2:end))];
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
    
    save(fullfile(savepath,'expe_merged_predicted_effects_accuracy.mat'),...
        'ndata','modnam','models','nmodel','xavg_fit','xstd_fit',...
        'parnam','npar','stab_param_nam','stab_param_vec',...
        'pcor_dat_avg','pcor_dat_std',...
        'pcor_mod_avg','pcor_mod_std',...
        'pcor_sim_avg','pswi_sim_avg',...
        'lastbef_sim_avg','firstaf_sim_avg')
end
%
fprintf('Optimal performance achieved: \n')
pbar = .9; HGT = 4;
for imodel = 1:nmodel
    % interpolate to find maximum
    xdat = stab_param_vec(imodel,:)';
    ydat = mean(pcor_sim_avg(:,:,imodel),1)';
    par_opt     = fminbnd(@(x)-interp1(xdat,ydat,x,'spline'),0,max(xdat));
    pcor_opt    = interp1(xdat,ydat,par_opt,'spline');
    pcor0       = interp1(xdat,ydat,0,'spline');
    
    figname = sprintf('Figure_8a_accuracy_simulation_model_%s',modnam{imodel});
    figure('Color','white','Name',figname)
    hold on
    xlim([min(xdat) max(xdat)]);
    ylim([0.5 0.9]);
    plot(par_opt*[1,1],ylim,'k-');
    plot(xdat,ydat,'-','LineWidth',1.5,'Color',rgb(imodel,:));
    plot(par_opt,pcor_opt,'ko','MarkerSize',10,'MarkerFaceColor',rgb(imodel,:));
    
    gain = (pcor_opt-pcor0)*100;
    fprintf('\t%s:\t%4.2f%% (= p0%+4.3f) for %s \t= %4.3f\n',models{imodel},pcor_opt*100,gain,stab_param_nam{imodel},par_opt)
    
    switch stab_param_nam{imodel}
        case 'lambda'
            set(gca,'XTick',0:.5:1.5,'XTickLabel',{'0','0.5','1.0','1.5'});
        case {'delta','beta'}
            set(gca,'XTick',0:1:3,'XTickLabel',{'0','1.0','2.0','3.0'});
        case {'plaps','epsi'}
            set(gca,'XTick',0:.5:1);
    end
    plot(xlim,pcor0*[1 1],'k:')
    xlabel(sprintf('%s',stab_param_nam{imodel}));
    set(gca,'YTick',0.5:0.1:0.9);
    ylabel('fraction correct');
    helper_ILL(gcf,pbar,HGT);
end
fprintf('\n')

% Observed and predicted effects of best-fitting model parameters on accuracy
imodel = 3; % conditional inference model
nsmp = 1e4; % number of samples for estimating peak value and s.e.
for regtyp = {'data','model'} % regress human data or model simulations
    switch regtyp{:}
        case 'data'
            mtyp = 'o';
            ydat = pcor_dat_avg;
            ystd = pcor_dat_std;
        case 'model'
            mtyp = 'd';
            ydat = pcor_mod_avg(:,imodel);
            ystd = pcor_mod_std(:,imodel);
    end
    xdat = cell2mat(xavg_fit(:,imodel));
    xstd = cell2mat(xstd_fit(:,imodel));
    % fit regression model
    modfun = 'y ~ b1 + b2*x1 + b3*x1^2 + b4*x2 + b5*x3 + b6*x3^2';
    b0 = cat(2,mean(ydat),zeros(1,5));
    opt = statset('nlinfit');
    opt.RobustWgtFun = 'bisquare';
    nlm = fitnlm(xdat,ydat,modfun,b0,'Options',opt);
    w = nlm.Robust.Weights;
    fprintf('<strong>Regression model for fraction correct (%s):</strong>\n',regtyp{:});
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
        switch regtyp{:}
            case 'data'
                figname = sprintf('Figure_8b_accuracy_multiple_regression_%s_%s_%s',regtyp{:},modnam{imodel},parnam{imodel}{ipar});
            case 'model'
                figname = sprintf('Figure_8c_accuracy_multiple_regression_%s_%s_%s',regtyp{:},modnam{imodel},parnam{imodel}{ipar});
        end
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
        switch regtyp{:}
            case 'data'
                figname = sprintf('Figure_8b_accuracy_multiple_regression_%s_%s_%s_coeff',regtyp{:},modnam{imodel},parnam{imodel}{ipar});
            case 'model'
                figname = sprintf('Figure_8c_accuracy_multiple_regression_%s_%s_%s_coeff',regtyp{:},modnam{imodel},parnam{imodel}{ipar});
        end
        figure('Color','white','Name',figname)
        pbar = 0.5; HGT = 1.5;
        hold on
        switch parnam{imodel}{ipar}
            case 'h'
                icoeff = 3;
            case 'siginf'
                icoeff = 4; ylim([-0.05,0]);
            case {'delta' 'plaps' 'lambda' 'epsi' 'beta'}
                icoeff = 6; ylim([-0.04,0]);
        end
        bavg = nlm.Coefficients.Estimate(icoeff);
        berr = nlm.Coefficients.SE(icoeff);
        xlim([0.6,1.4]);
        
        bar(1,bavg,0.5,'EdgeColor',rgb(imodel,:),'FaceColor',0.5*(rgb(imodel,:)+1),'LineWidth',1);
        plot([1,1],bavg+berr*[-1,+1],'k-');
        set(gca,'XAxisLocation','top','XTick',[]);
        switch parnam{imodel}{ipar}
            case {'siginf', 'delta'}
                set(gca,'YTick',-0.04:0.02:0,'YTickLabel',{'-0.04' '-0.02','0'});
        end
        xlabel(sprintf('%s',nlm.CoefficientNames{icoeff-1}));
        helper_ILL(gcf,pbar,HGT,false);
    end
    fprintf('\n')
end