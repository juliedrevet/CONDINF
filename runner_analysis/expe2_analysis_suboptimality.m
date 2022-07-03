%% EXPERIMENT 2 RUNNER ANALYSIS SUBOPTIMALITY
% clear workspace
clearvars

savepath = '../results/';
datapath = '../data_table';

addpath '../helper_plot';
addpath '../matlab_functions';
%% CHARACTERIZING THE SUBOPTIMALITY OF STATISTICAL INFERENCE - REPLICATION (Figure 5d-f)
% reversal behaviors, overall accuracy, overall switch rate
if exist(fullfile(savepath,'expe2_suboptimality.mat'),'file')
    load(fullfile(savepath,'expe2_suboptimality.mat'))
    ndata = size(pcor_dat,1);
else
    fprintf('<strong>Compute reversal behaviors, overall accuracy, overall switch rate</strong>\n')
    % load data table
    load('./data_table/expe2_data.mat','expe2_data');
    dat = expe2_data;
    clear('expe2_data');
    % get list of participants
    participants = unique(dat.participant);
    ndata = numel(participants);
    
    % PARTICIPANTS
    fprintf('\tparticipants...')
    pcor_dat = nan(ndata,1);
    prep_dat = nan(ndata,1);
    crev_dat = nan(ndata,8);
    crep_dat = nan(ndata,8);
    crev_w_dat = nan(ndata,8,3);
    crep_w_dat = nan(ndata,8,3);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);   % trial number in block
        s = dat.stim(ifilt);    % stimulus (+1/-1)
        w = dat.wght(ifilt);    % stimulus weight or strength (1/2/3)
        r = dat.resp(ifilt);    % participant's response (+/-1)
        
        pcor_dat(idata) = nanmean(s==r);
        rep = [nan; r(1:end-1) == r(2:end)];
        rep(t==1) = nan;
        prep_dat(idata) = nanmean(rep);
        c = getc(t,s,r,w);
        crev_dat(idata,:) = c.rev_avg;
        crep_dat(idata,:) = c.rep_avg;
        crev_w_dat(idata,:,:) = c.rev_w_avg; % one per stimulus strength on the trial after reversal
        crep_w_dat(idata,:,:) = c.rep_w_avg;
    end
    fprintf('done.\n')
    
    % OPTIMAL MODEL
    fprintf('\toptimal model...')
    pcor_opt = nan(ndata,1);
    prep_opt = nan(ndata,1);
    crev_opt = nan(ndata,8);
    crep_opt = nan(ndata,8);
    crev_w_opt = nan(ndata,8,3);
    crep_w_opt = nan(ndata,8,3);
    h_opt = nan(ndata,1);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);           % trial number in block
        s = dat.stim(ifilt);            % stimulus (+1/-1)
        w = dat.wght(ifilt);            % stimulus weight or strength (1/2/3)
        r = dat.resp(ifilt);            % participant's response (+/-1)
        isrev = dat.reversal(ifilt);    % reversal trial
        h_true = mean(isrev);           % true hazard rate
        
        % configure model simulation
        cfg_sim        = [];
        cfg_sim.t      = t;
        cfg_sim.s      = s;
        cfg_sim.w      = w;
        cfg_sim.modtype  = 'flatinf';
        cfg_sim.pfal   = [0.30 0.20 0.10];
        cfg_sim.h      = h_true;    % true hazard rate
        cfg_sim.siginf = 0;         % no inference noise
        cfg_sim.sigsel = 0;         % no selection noise
        cfg_sim.nsmp   = 1e4;       % number of simulations
        
        out_sim = sim_model_expe2(cfg_sim);
        
        pcor_opt(idata) = out_sim.pcor;
        prep_opt(idata) = out_sim.prep;
        crev_opt(idata,:) = out_sim.c.rev_avg;
        crep_opt(idata,:) = out_sim.c.rep_avg;
        crev_w_opt(idata,:,:) = out_sim.c.rev_w_avg;
        crep_w_opt(idata,:,:) = out_sim.c.rep_w_avg;
        clearvars out_sim;
    end
    fprintf('done.\n');
    
    % SUBOPTIMAL NOISY INFERENCE MODEL
    fprintf('\tsuboptimal noisy inference model...')
    load(fullfile(savepath,'expe2_results_fit_noisy_inference_model.mat'),'out_fit')
    % extract best-fitting parameters
    h = cellfun(@(s)getfield(s,'h'),out_fit);
    siginf = cellfun(@(s)getfield(s,'siginf'),out_fit);
    pcor_sub = nan(ndata,1);
    prep_sub = nan(ndata,1);
    crev_sub = nan(ndata,8);
    crep_sub = nan(ndata,8);
    crev_w_sub = nan(ndata,8,3);
    crep_w_sub = nan(ndata,8,3);
    for idata = 1:ndata
        fprintf('.')
        % filter participant
        ifilt = dat.participant == participants(idata);
        % get data
        t = dat.trial(ifilt);	% trial number in block
        s = dat.stim(ifilt);	% stimulus (+1/-1)
        w = dat.wght(ifilt);	% stimulus weight or strength (1/2/3)
        r = dat.resp(ifilt);	% participant's response (+/-1)
        
        % configure model simulation
        cfg_sim        = [];
        cfg_sim.t      = t;
        cfg_sim.s      = s;
        cfg_sim.w      = w;
        cfg_sim.modtype  = 'nonstab';
        cfg_sim.pfal   = [0.30 0.20 0.10];  % false-alarm rate
        cfg_sim.h      = h(idata);          % perceived hazard rate
        cfg_sim.siginf = siginf(idata);     % fitted inference noise
        cfg_sim.sigsel = 0;                 % no selection noise
        cfg_sim.nsmp   = 1e4;               % number of simulations
        
        out_sim = sim_model_expe2(cfg_sim);
        
        pcor_sub(idata) = out_sim.pcor;
        prep_sub(idata) = out_sim.prep;
        crev_sub(idata,:) = out_sim.c.rev_avg;
        crep_sub(idata,:) = out_sim.c.rep_avg;
        crev_w_sub(idata,:,:) = out_sim.c.rev_w_avg;
        crep_w_sub(idata,:,:) = out_sim.c.rep_w_avg;
        clearvars out_sim;
    end
    fprintf('done.\n');
    save(fullfile(savepath,'expe2_suboptimality.mat',...
        'pcor_dat','prep_dat','crev_dat','crep_dat','crev_w_dat','crep_w_dat',...
        'pcor_opt','prep_opt','crev_opt','crep_opt','crev_w_opt','crep_w_opt','h_true',...
        'pcor_sub','prep_sub','crev_sub','crep_sub','crev_w_sub','crep_w_sub','h'))
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

[H,P,~,STATS] = ttest(pcor_opt,pcor_dat);
fprintf('\n<strong>overall accuracy optimal inference versus true data</strong>\n\t H = %d ; t = %5.3f ; p > %5.3f\n',H,STATS.tstat, floor(1000*P)/1000);
[H,P,~,STATS] = ttest(pcor_sub,pcor_dat);
fprintf('\n<strong>overall accuracy suboptimal noisy inference versus true data</strong>\n\t H = %d ; t = %5.3f ; p < %5.3f\n',H,STATS.tstat, ceil(1000*P)/1000);

[H,P,~,STATS] = ttest(1-prep_opt,1-prep_dat);
fprintf('\n<strong>overall switch rate optimal inference versus true data</strong>\n\t H = %d ; t = %5.3f ; p > %5.3f\n',H,STATS.tstat, floor(1000*P)/1000);
[H,P,~,STATS] = ttest(1-prep_sub,1-prep_dat);
fprintf('\n<strong>overall switch rate suboptimal noisy inference versus true data</strong>\n\t H = %d ; t = %5.3f ; p < %5.3f\n',H,STATS.tstat, ceil(1000*P)/1000);

[H,P,~,STATS] = ttest(h-h_true);
fprintf('\n<strong>Suboptimal inference</strong>\n')
fprintf('\thazard rate h\t\t %5.3f +/- %5.3f\n', mean(h),std(h)/sqrt(ndata));
fprintf('\ttrue hazard rate h\t %5.3f\n',h_true);
fprintf('\tperceived vs true hazard rate\t H = %d ; t = %5.3f ; p > %5.3f\n',H,STATS.tstat, floor(1000*P)/1000);

%% REVERSAL BEHAVIOR (Figure 5d-e)
rgb_dat0 = [0.5 0.5 0.5];
rgb_opt0 = [157 44 41]/255;
rgb_sub0 = [190 215 231]/255;

rgb_dat = nan(3,3); rgb_opt = nan(3,3); rgb_sub = nan(3,3);
for istr = 1:3
    rgb_dat(istr,:) = (rgb_dat0+(3-istr))/(4-istr);
    rgb_opt(istr,:) = (rgb_opt0+(3-istr))/(4-istr);
    rgb_sub(istr,:) = (rgb_sub0+(3-istr))/(4-istr);
end

% plot optimal model versus human data
cfg = [];
cfg.figname = 'Figure_5d_expe2_REV_optimal';
cfg.rgb = rgb_opt;
cfg.rgb_dat = rgb_dat;
cfg.yavg = squeeze(mean(crev_w_opt,1));
cfg.Yavg = squeeze(mean(crev_w_dat,1));
cfg.Yerr = squeeze(std(crev_w_dat,[],1))/sqrt(ndata);
helper_plot_rev(cfg)

cfg = [];
cfg.figname = 'Figure_5d_expe2_SWI_optimal';
cfg.rgb = rgb_opt;
cfg.rgb_dat = rgb_dat;
cfg.yavg = 1-squeeze(mean(crep_w_opt,1));
cfg.Yavg = 1-squeeze(mean(crep_w_dat,1));
cfg.Yerr = squeeze(std(crep_w_dat,[],1))/sqrt(ndata);
cfg.ylim = [0 0.6]; 
helper_plot_rep(cfg)

% plot noisy inference model versus human data
cfg = [];
cfg.figname = 'Figure_5e_expe2_REV_noisy_inference';
cfg.rgb = rgb_sub;
cfg.rgb_dat = rgb_dat;
cfg.yavg = squeeze(mean(crev_w_sub,1));
cfg.yerr = squeeze(std(crev_w_sub,[],1))/sqrt(ndata);
cfg.Yavg = squeeze(mean(crev_w_dat,1));
cfg.Yerr = squeeze(std(crev_w_dat,[],1))/sqrt(ndata);
helper_plot_rev(cfg)

cfg = [];
cfg.figname = 'Figure_5e_expe2_SWI_noisy_inference';
cfg.rgb = rgb_sub;
cfg.rgb_dat = rgb_dat;
cfg.yavg = 1-squeeze(mean(crep_w_sub,1));
cfg.yerr = squeeze(std(crep_w_sub,[],1))/sqrt(ndata);
cfg.Yavg = 1-squeeze(mean(crep_w_dat,1));
cfg.Yerr = squeeze(std(crep_w_dat,[],1))/sqrt(ndata);
cfg.ylim = [0 0.6]; 
helper_plot_rep(cfg)

%% OVERALL BEHAVIOR (Figure 5f)
% compare overall accuracy
cfg = [];
cfg.figname = 'Figure_5f_expe2_accuracy';
cfg.rgb = [rgb_opt0;rgb_sub0];
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
cfg.figname = 'Figure_5f_expe2_switch_rate';
cfg.rgb = [rgb_opt0;rgb_sub0];
cfg.x = 1-[prep_opt prep_sub];
cfg.xdat = 1-prep_dat;
cfg.ylim = [0 .5];
cfg.ylabel = 'overall switch rate';
cfg.xlabel = 'inference';
cfg.xticklabel = {'optimal' 'noisy'};
cfg.HGT = 4; cfg.pbar = .9;
helper_plot_overall_bars(cfg);

