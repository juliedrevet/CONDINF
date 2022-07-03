%% EXPERIMENT 2 RUNNER ANALYSIS STABILIZATION
% needed on the path for running analysis:
%       * simple_mixed_anova.m (Calpette, L. (2022) https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)

% clear workspace
clearvars
addpath '../helper_plot';
addpath '../matlab_functions';

savepath = '../results/';
rgb = [ ...
    238,197,141; ... % senbias
    048,102,158; ... % inflaps
    134,198,226; ... % infdisc
    171,131,171; ... % selbias
    099,074,135  ... % selepsi
    ]/255;

%% VALIDATION PREDICTIONS CONDITIONAL INFERENCE (Figure 6a)
fprintf('<strong>Validation of specific predictions of conditional inference</strong>\n')
load(fullfile(savepath,'expe2_results_fit_conditional_inference_independent_thresholds.mat'),'out_fit')
ndata = size(out_fit,1);

rgb_mod = rgb(3,:); % conditional inference model color
xavg_fit = cell2mat(cellfun(@(s)getfield(s,'xavg'),out_fit,'UniformOutput',0));
delta_vec = xavg_fit(:,3:5);
fprintf('ANOVA\n')

tbl = simple_mixed_anova(delta_vec,[],{'reliability_threshold'}); disp(tbl);
fprintf('eta_squared = %.5f\n\n',tbl.SumSq(3)/sum(tbl.SumSq(3:4)))

% compute discard rate from reliability threshold
ndata = size(delta_vec,1);
nstr = size(delta_vec,2);

pfal = [0.30 0.20 0.10]';       % false-alarm rate
sigsen = 1/norminv(1-pfal(2));  % sensory noise
strsen = norminv(1-pfal)/norminv(1-pfal(2)); % stimulus sensory strengths

ph0_vec = nan(ndata,nstr);  % probability of hits
pf0_vec = nan(ndata,nstr);  % probability of false-alarms
pl0_vec = nan(ndata,nstr);  % probability of discard

for istr = 1:nstr
    ph0_vec(:,istr) = 1-normcdf(+delta_vec(:,istr),strsen(istr)*2/sigsen^2,2/sigsen);
    pf0_vec(:,istr) = normcdf(-delta_vec(:,istr),strsen(istr)*2/sigsen^2,2/sigsen);
    pl0_vec(:,istr) = 1-ph0_vec(:,istr)-pf0_vec(:,istr);
end
fprintf('ANOVA\n')

tbl = simple_mixed_anova(pl0_vec,[],{'discard_rate'}); disp(tbl);
fprintf('eta_squared = %.5f\n\n',tbl.SumSq(3)/sum(tbl.SumSq(3:4)))

fprintf('best-fitting independent reliability thresholds\n')
fprintf('\t* delta1 = %.3f +/- %.3f\n\t* delta2 = %.3f +/- %.3f\n\t* delta3 = %.3f +/- %.3f\n',...
        mean(delta_vec(:,1)),std(delta_vec(:,1))/sqrt(ndata),...
        mean(delta_vec(:,2)),std(delta_vec(:,2))/sqrt(ndata),...
        mean(delta_vec(:,3)),std(delta_vec(:,3))/sqrt(ndata))

fprintf('\ncorresponding discard rates\n')
fprintf('\t* pdisc(70%%) = %.3f +/- %.3f\n\t* pdisc(80%%) = %.3f +/- %.3f\n\t* pdisc(90%%) = %.3f +/- %.3f\n',...
        mean(pl0_vec(:,1)),std(pl0_vec(:,1))/sqrt(ndata),...
        mean(pl0_vec(:,2)),std(pl0_vec(:,2))/sqrt(ndata),...
        mean(pl0_vec(:,3)),std(pl0_vec(:,3))/sqrt(ndata))
    
% violinplot of independently fitted reliability thresholds   
figname = 'Figure_6a_independant_reliability_thresholds';
figure('Color','white','Name',figname);
pbar = 1.2; HGT = 4;
y = delta_vec;
hold on
nstim = size(y,2);
xlim([0.4,nstim+.6]); ylim([0,3]);
for istim = 1:nstim
    x = y(:,istim);
    pos = istim;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    [pk,~,u] = ksdensity(x,xk);
    while true
        z = findpeaks(pk,xk);
        if numel(z) <= 1
            break
        else
            u = u*1.1;
            pk = ksdensity(x,xk,'BandWidth',u);
        end
    end
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    xmed = interp1(xk,pk,median(x));
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb_mod+1),'EdgeColor','none','FaceAlpha',0.2+istim/5);
    plot(pos+pk*wid,xk,'-','Color',rgb_mod,'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb_mod,'LineWidth',0.5);
    plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-');
    for i = 1:numel(x)
        plot(pos+jit(i)*str(i)*wid,x(i),'wo','MarkerFaceColor',rgb_mod,'MarkerSize',4,'LineWidth',0.5);
    end
    plot(pos,mean(x),'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',0.5);
    plot(pos+xmed*wid*[-1;+1],median(x).*[1 1],'Color',rgb_mod,'LineWidth',0.5);
end
set(gca,'XTick',1:nstim,'XTickLabel',{'70%','80%','90%'},'YTick',0:0.5:3);
xlabel('stimulus strength')
ylabel('reliability threshold \delta')
helper_ILL(gcf,pbar,HGT);

% violinplot of corresponding discard rates  
figname = 'Figure_6a_independant_discard_rates';
figure('Color','white','Name',figname)
y = pl0_vec;
hold on
nstim = size(y,2);
xlim([0.4,nstim+.6]); ylim([0,1]);
for istim = 1:nstim
    x = y(:,istim);
    pos = istim;
    wid = 0.4;
    xk = linspace(min(x),max(x),100);
    [pk,~,u] = ksdensity(x,xk);
    while true
        z = findpeaks(pk,xk);
        if numel(z) <= 1
            break
        else
            u = u*1.1;
            pk = ksdensity(x,xk,'BandWidth',u);
        end
    end
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    xmed = interp1(xk,pk,median(x));
    jit = linspace(-1,+1,numel(x));
    patch([pos+pk*wid,fliplr(pos-pk*wid)],[xk,fliplr(xk)],0.5*(rgb_mod+1),'EdgeColor','none','FaceAlpha',0.2+istim/5);
    plot(pos+pk*wid,xk,'-','Color',rgb_mod,'LineWidth',0.5);
    plot(pos-pk*wid,xk,'-','Color',rgb_mod,'LineWidth',0.5);
    plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-');
    for i = 1:numel(x)
        plot(pos+jit(i)*str(i)*wid,x(i),'wo','MarkerFaceColor',rgb_mod,'MarkerSize',4,'LineWidth',0.5);
    end
    plot(pos,mean(x),'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',0.5);
    plot(pos+xmed*wid*[-1;+1],median(x).*[1 1],'Color',rgb_mod,'LineWidth',0.5);
end
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'XTick',1:nstim,'XTickLabel',{'70%','80%','90%'},'YTick',0:0.2:1);
xlabel('stimulus strength')
ylabel('discard rate')
helper_ILL(gcf,pbar,HGT);

% RELATION WITH CONFIDENCE
load(fullfile(savepath,'expe2_results_fit_confidence.mat'));
% get best-fitting confidence threshold parameter:
xavg_fit = cell2mat(cellfun(@(s) s.xavg ,out_fit,'UniformOutput',false));
xstd_fit = cell2mat(cellfun(@(s) s.xstd ,out_fit,'UniformOutput',false));
thrcon_avg = xavg_fit(:,3);
thrcon_std = xstd_fit(:,3);
% get best-fitting reliability threshold parameter:
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');
imodel = 3; % conditional inference model
xavg_fit = cell2mat(cellfun(@(s) s.xavg ,out_fit(:,imodel),'UniformOutput',false));
xstd_fit = cell2mat(cellfun(@(s) s.xstd ,out_fit(:,imodel),'UniformOutput',false));
delta_avg = xavg_fit(:,3);
delta_std = xstd_fit(:,3);

% plot correlation between confidence threshold and reliability threshold
fprintf('\nRelation with confidence\n');
cfg = [];
cfg.figname = 'Figure_6b_correlation_confidence_reliability';
cfg.xavg = thrcon_avg;
cfg.xstd = thrcon_std;
cfg.yavg = delta_avg;
cfg.ystd = delta_std;
cfg.xlim = [-3,5]; cfg.xtick = -2:2:+4; cfg.xticklabel = {'-2.0','0','2.0','4.0'};
cfg.xlabel = sprintf('confidence threshold delta_{c}');
cfg.ylim = [0 2.5]; cfg.ytick = 0:1:2.5; cfg.yticklabel = {'0','1.0','2.0'};
cfg.ylabel = sprintf('reliability threshold delta');
cfg.pbar = 1; cfg.HGT = 4;
cfg.rgb = rgb(imodel,:);
cfg.ci = true;
helper_plot_corr(cfg);

[~,~,RLO,RUP]=corrcoef(cfg.xavg ,cfg.yavg);
fprintf('\t\tCI = [%+5.3f %+5.3f]\n', RLO(1,2),RUP(1,2));

%% VARIANCE OF CONFIDENCE THRESHOLD EXPLAINED BY EACH MODEL PARAMETER (Figure 6b)
% get best-fitting conditional inference parameters:
load(fullfile(savepath,'expe2_results_fit_stabilizing_models.mat'),'out_fit');
imodel = 3; % conditional inference model
xavg_fit = cell2mat(cellfun(@(s) s.xavg ,out_fit(:,imodel),'UniformOutput',false));
xstd_fit = cell2mat(cellfun(@(s) s.xstd ,out_fit(:,imodel),'UniformOutput',false));
ndata = size(xavg_fit,1);
npar = size(xavg_fit,2);
parnam = out_fit{1,imodel}.xnam;

% compute bootstrap resamples
nres = 1e4; % number of bootstrap resamples
rho_res = nan(nres,3); % correlation strength
rsq_res = nan(nres,3); % variance explained
for ires = 1:nres
    idata = randsample(ndata,ndata,true);
    for ipar = 1:npar
        x = thrcon_avg(idata,:);
        y = xavg_fit(idata,ipar);
        rho_tmp = corr(x,y);
        rho_res(ires,ipar) = rho_tmp;
        rsq_res(ires,ipar) = rho_tmp^2;
    end
end
% compute bootstraped statistics
rsq = median(rsq_res,1); % median
rsq_qt = nan(3,2); % 1st/3rd quartiles
for ipar = 1:npar
    rsq_qt(ipar,:) = quantile(rsq_res(:,ipar),[1,3]/4);
end
% compute exceedance probabilities for variance explained
pexc = nan(3,1);
[~,imax] = max(rsq_res,[],2);
for ipar = 1:npar
    pexc(ipar) = mean(imax == ipar);
end

% plot variance of confidence threshold explained by each model parameter
cfg = [];
cfg.figname = 'Figure_6b_confidence_threshold_variance_explained';
[~,Imax] = max(pexc);
cfg.rgb = repmat((1+rgb(imodel,:))/2,npar,1);
cfg.rgb(Imax,:) = rgb(imodel,:);
cfg.x   = rsq;
cfg.xqt = rsq_qt;
cfg.ylim = [0 .3];
cfg.ylabel = 'variance explained';
cfg.xlabel = 'model parameter';
cfg.xticklabel = parnam;
cfg.HGT = 4; cfg.pbar = 1;
helper_plot_overall_bars(cfg);
fprintf('Parameter %s: \texceedance p = %.3f\n\n',parnam{Imax},pexc(Imax));
   