%% EXPERIMENT 1 RUNNER ANALYSIS SUBOPTIMALITY
% clear workspace
clearvars
savepath = '../results/';
datapath = '../data_table';

addpath '../helper_plot';
addpath '../matlab_functions';

%% CHARACTERIZING THE SUBOPTIMALITY OF STATISTICAL INFERENCE (Figure 1)
% reversal behaviors, overall accuracy, overall switch rate
if exist(fullfile(savepath,'expe1_suboptimality.mat'),'file')
    load(fullfile(savepath,'expe1_suboptimality.mat'))
    ndata = size(pcor_dat,1);
else
    fprintf('<strong>Compute reversal behaviors, overall accuracy, overall switch rate</strong>\n')
    % load data table
    load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
    dat = expe1_data;
    clear('expe1_data');
    % get list of participants
    participants = setdiff(unique(dat.participant),[33]);
    ndata = numel(participants);
    
    % PARTICIPANTS
    fprintf('\tparticipants...')
    pcor_dat = nan(ndata,1);
    prep_dat = nan(ndata,1);
    crev_dat = nan(ndata,8);
    crep_dat = nan(ndata,8);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);   % trial number in block
        s = dat.stim(ifilt);    % stimulus category (+/-1)
        r = dat.resp(ifilt);    % participant's response (+/-1)
        
        pcor_dat(idata) = nanmean(s==r);
        rep = [nan; r(1:end-1) == r(2:end)];
        rep(t==1) = nan;
        prep_dat(idata) = nanmean(rep);
        c = getc(t,s,r);
        crev_dat(idata,:) = c.rev_avg;
        crep_dat(idata,:) = c.rep_avg;
    end
    fprintf('done.\n');
    
    % OPTIMAL MODEL
    fprintf('\toptimal model...')
    pcor_opt = nan(ndata,1);
    prep_opt = nan(ndata,1);
    crev_opt = nan(ndata,8);
    crep_opt = nan(ndata,8);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);           % trial number in block
        s = dat.stim(ifilt);            % stimulus (+/-1)
        r = dat.resp(ifilt);            % participant's response (+/-1)
        isrev = dat.reversal(ifilt);    % reversal trial
        h_true = mean(isrev);           % true hazard rate
        
        % configure model simulation
        cfg_sim        = [];
        cfg_sim.t      = t;
        cfg_sim.s      = s;
        cfg_sim.modtype  = 'flatinf';
        cfg_sim.pfal   = 0.20;      % false-alarm rate
        cfg_sim.h      = h_true;    % true hazard rate
        cfg_sim.siginf = 0;         % no inference noise
        cfg_sim.sigsel = 0;         % no selection noise
        cfg_sim.nsmp   = 1e4;       % number of simulations
        
        out_sim = sim_model(cfg_sim);
        
        pcor_opt(idata) =   out_sim.pcor;
        prep_opt(idata) =   out_sim.prep;
        crev_opt(idata,:) = out_sim.c.rev_avg;
        crep_opt(idata,:) = out_sim.c.rep_avg;
        clearvars out_sim;
    end
    fprintf('done.\n');
    
    % SUBOPTIMAL NOISY INFERENCE MODEL
    fprintf('\tsuboptimal noisy inference model...')
    load(fullfile(savepath,'expe1_results_fit_noisy_inference_model.mat'))
    % extract best-fitting parameters
    h = cellfun(@(s)getfield(s,'h'),out_fit);
    siginf = cellfun(@(s)getfield(s,'siginf'),out_fit);
    pcor_sub = nan(ndata,1);
    prep_sub = nan(ndata,1);
    crev_sub = nan(ndata,8);
    crep_sub = nan(ndata,8);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);	% trial number in block
        s = dat.stim(ifilt);    % stimulus (+/-1)
        r = dat.resp(ifilt);    % participant's response (+/-1)
        
        % configure model simulation
        cfg_sim        = [];
        cfg_sim.t      = t;
        cfg_sim.s      = s;
        cfg_sim.modtype  = 'nonstab';
        cfg_sim.pfal   = 0.20;              % false-alarm rate
        cfg_sim.h      = h(idata);          % perceived hazard rate
        cfg_sim.siginf = siginf(idata);     % fitted inference noise
        cfg_sim.sigsel = 0;                 % no selection noise
        cfg_sim.nsmp   = 1e4;               % number of simulations
        
        out_sim = sim_model(cfg_sim);
        
        pcor_sub(idata) = out_sim.pcor;
        prep_sub(idata) = out_sim.prep;
        crev_sub(idata,:) = out_sim.c.rev_avg;
        crep_sub(idata,:) = out_sim.c.rep_avg;
        clearvars out_sim;
    end
    fprintf('done.\n');
    save(fullfile(savepath,'expe1_suboptimality.mat'),...
        'pcor_dat','prep_dat','crev_dat','crep_dat',...
        'pcor_opt','prep_opt','crev_opt','crep_opt','h_true',...
        'pcor_sub','prep_sub','crev_sub','crep_sub','h')
end

fprintf('\n<strong>Participants</strong>\n')
fprintf('\toverall accuracy\t %5.3f +/- %5.3f\n', mean(pcor_dat),std(pcor_dat)/sqrt(ndata));
fprintf('\toverall switch rate\t %5.3f +/- %5.3f\n', 1-mean(prep_dat),std(prep_dat)/sqrt(ndata));

fprintf('\n<strong>Bayes-optimal inference</strong>\n')
fprintf('\toverall accuracy\t %5.3f +/- %5.3f\n', mean(pcor_opt),std(pcor_opt)/sqrt(ndata));
fprintf('\toverall switch rate\t %5.3f +/- %5.3f\n', 1-mean(prep_opt),std(prep_opt)/sqrt(ndata));

fprintf('\n<strong>Suboptimal noisy inference</strong>\n')
fprintf('\toverall accuracy\t %5.3f +/- %5.3f\n', mean(pcor_sub),std(pcor_sub)/sqrt(ndata));
fprintf('\toverall switch rate\t %5.3f +/- %5.3f\n', 1-mean(prep_sub),std(prep_sub)/sqrt(ndata));

[H,P,CI,STATS] = ttest(pcor_opt,pcor_dat); effd = STATS.tstat/sqrt(size(pcor_dat,1));
if H == 1, psign = '<';else,psign = '>';end
fprintf('\n<strong>overall accuracy optimal inference versus true data</strong>\n\t ')
fprintf('H = %d ; t = %5.3f ; p %s %5.3f ; Cohen''s d = %5.3f ; CI = [%+5.3f %+5.3f]\n',H,STATS.tstat,psign,floor(1000*P)/1000,effd,CI(1),CI(2));

[H,P,CI,STATS] = ttest(pcor_sub,pcor_dat); effd = STATS.tstat/sqrt(size(pcor_dat,1));
if H == 1, psign = '<';else,psign = '>';end
fprintf('\n<strong>overall accuracy suboptimal noisy inference versus true data</strong>\n\t ')
fprintf('H = %d ; t = %5.3f ; p %s %5.3f ; Cohen''s d = %5.3f ; CI = [%+5.3f %+5.3f]\n',H,STATS.tstat,psign,ceil(1000*P)/1000,effd,CI(1),CI(2));

[H,P,CI,STATS] = ttest(1-prep_opt,1-prep_dat); effd = STATS.tstat/sqrt(size(1-prep_dat,1));
if H == 1, psign = '<';else,psign = '>';end
fprintf('\n<strong>overall switch rate optimal inference versus true data</strong>\n\t ')
fprintf('H = %d ; t = %5.3f ; p %s %5.3f ; Cohen''s d = %5.3f ; CI = [%+5.3f %+5.3f]\n',H,STATS.tstat,psign,floor(1000*P)/1000,effd,CI(1),CI(2));

[H,P,CI,STATS] = ttest(1-prep_sub,1-prep_dat); effd = STATS.tstat/sqrt(size(1-prep_dat,1));
if H == 1, psign = '<';else,psign = '>';end
fprintf('\n<strong>overall switch rate suboptimal noisy inference versus true data</strong>\n\t ')
fprintf('H = %d ; t = %5.3f ; p %s %5.3f ; Cohen''s d = %5.3f ; CI = [%+5.3f %+5.3f]\n',H,STATS.tstat,psign,ceil(1000*P)/1000,effd,CI(1),CI(2));

[H,P,CI,STATS] = ttest(h-h_true); effd = STATS.tstat/sqrt(size(h,1));
if H == 1, psign = '<';else,psign = '>';end
fprintf('\n<strong>Suboptimal inference</strong>\n')
fprintf('\thazard rate h\t\t %5.3f +/- %5.3f\n', mean(h),std(h)/sqrt(ndata));
fprintf('\ttrue hazard rate h\t %5.3f\n',h_true);
fprintf('\tperceived vs true hazard rate\t ')
fprintf('H = %d ; t = %5.3f ; p %s %5.3f ; Cohen''s d = %5.3f ; CI = [%+5.3f %+5.3f]\n\n',H,STATS.tstat,psign,floor(1000*P)/1000,effd,CI(1),CI(2));

%% REVERSAL BEHAVIOR (Figure 1e - Figure 1f)
rgb_dat = [0.5 0.5 0.5];
rgb_opt = [157 44 41]/255;
rgb_sub = [190 215 231]/255;

% plot optimal model versus human data
cfg = [];
cfg.figname = 'Figure_1e_expe1_REV_optimal';
cfg.rgb = rgb_opt;
cfg.rgb_dat = rgb_dat;
cfg.yavg = mean(crev_opt,1)';
cfg.Yavg = mean(crev_dat,1)';
cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
helper_plot_rev(cfg)

cfg = [];
cfg.figname = 'Figure_1e_expe1_SWI_optimal';
cfg.rgb = rgb_opt;
cfg.rgb_dat = rgb_dat;
cfg.yavg = 1-mean(crep_opt,1)';
cfg.Yavg = 1-mean(crep_dat,1)';
cfg.Yerr = std(crep_dat,[],1)'/sqrt(ndata);
helper_plot_rep(cfg)

% plot noisy inference model versus human data
cfg = [];
cfg.figname = 'Figure_1f_expe1_REV_noisy_inference';
cfg.rgb = rgb_sub;
cfg.rgb_dat = rgb_dat;
cfg.yavg = mean(crev_sub,1)';
cfg.yerr = std(crev_sub,[],1)'/sqrt(ndata);
cfg.Yavg = mean(crev_dat,1)';
cfg.Yerr = std(crev_dat,[],1)'/sqrt(ndata);
helper_plot_rev(cfg)

cfg = [];
cfg.figname = 'Figure_1f_expe1_SWI_noisy_inference';
cfg.rgb = rgb_sub;
cfg.rgb_dat = rgb_dat;
cfg.yavg = 1-mean(crep_sub,1)';
cfg.yerr = std(crep_sub,[],1)'/sqrt(ndata);
cfg.Yavg = 1-mean(crep_dat,1)';
cfg.Yerr = std(crep_dat,[],1)'/sqrt(ndata);
helper_plot_rep(cfg)

%% OVERALL BEHAVIOR (Figure 1g)
% compare overall accuracy
cfg = [];
cfg.figname = 'Figure_1g_expe1_accuracy';
cfg.rgb = [rgb_opt;rgb_sub];
cfg.x = [pcor_opt pcor_sub];
cfg.xdat = pcor_dat;
cfg.ylim = [.5 1];
cfg.ylabel = 'overall accuracy';
cfg.xlabel = 'inference';
cfg.xticklabel = {'optimal' 'noisy'};
cfg.HGT = 4; cfg.pbar = .9;
helper_plot_overall_bars(cfg);

% compare overall switch rate
cfg = [];
cfg.figname = 'Figure_1g_expe1_switch_rate';
cfg.rgb = [rgb_opt;rgb_sub];
cfg.x = 1-[prep_opt prep_sub];
cfg.xdat = 1-prep_dat;
cfg.ylim = [0 .5];
cfg.ylabel = 'overall switch rate';
cfg.xlabel = 'inference';
cfg.xticklabel = {'optimal' 'noisy'};
cfg.HGT = 4; cfg.pbar = .9;
helper_plot_overall_bars(cfg);
