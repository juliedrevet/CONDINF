%% RUNNER ANALYSIS EXPERIMENT 2
% needed on the path for running analysis:
%       * simple_mixed_anova.m (Calpette, L. (2022) https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)

addpath '../helper_plot';
addpath '../matlab_functions';

fprintf('<strong>***\tEXPERIMENT 2\t***</strong>\n\n')

%% Relation between decision confidence and accuracy in titration trials (Figure 5c)
clearvars;

datapath = '../data_table';
fprintf('<strong>Relation between decision confidence and accuracy in titration trials</strong>\n')
% load data
load(fullfile(datapath,'expe2_titration.mat'),'expe2_titration');
dat = expe2_titration;
clear('expe2_titration');

% get list of participants
participants = unique(dat.participant);
ndata = numel(participants);

dat.bins = nan(size(dat,1),1);

meanconf = nan(ndata,7);
meancorr = nan(ndata,7);
meanconf_cor = nan(ndata,7);
meanconf_err = nan(ndata,7);
meancorr_conf = nan(ndata,7);
meancorr_uncf = nan(ndata,7);

fprintf('Computing fraction correct and relative confidence depending on binned marble luminance\n\t...')
for idata = 1:ndata
    fprintf('.')
    % filter participant
    ifilt = dat.participant == participants(idata);
    
    % get data
    blck = dat.block(ifilt); % titration block
    prop = dat.prop(ifilt);  % marble proportion ([0-100%] w.r.t. light)
    conf = dat.conf(ifilt);  % confidence (+1/-1)
    resp = dat.resp(ifilt);  % response (+1/-1)
    corr = dat.resp(ifilt) == dat.staircase(ifilt);
    bins = dat.bins(ifilt);
    
    conf(isnan(conf)) = 0;
    
    ifilt = dat.participant == participants(idata) & dat.trial == 1;
    xmax = dat.xmax(ifilt,:);
    xmin = dat.xmin(ifilt,:);
    
    xbins = [fliplr(xmin) xmax];
    
    for iblck = 1:16
        bins(blck == iblck) = discretize(prop(blck == iblck) ,[0 xbins(iblck,:) 1]);
    end
    
    % performance depending on confidence
    for ibin = 1:7
        meancorr(idata,ibin) = nanmean(corr(bins==ibin));
        meancorr_conf(idata,ibin) = nanmean(corr(bins==ibin & conf>0));
        meancorr_uncf(idata,ibin) = nanmean(corr(bins==ibin & conf<0));
    end
    
    % confidence depending on performance
    conf = conf-nanmean(conf);
    for ibin = 1:7
        meanconf(idata,ibin) = nanmean(conf(bins==ibin));
        meanconf_cor(idata,ibin) = nanmean(conf(bins==ibin & corr));
        meanconf_err(idata,ibin) = nanmean(conf(bins==ibin & ~corr));
    end
end
fprintf('done.\n')

% STATS
fprintf('ANOVA\n')

tbl = simple_mixed_anova(cat(3,meancorr_conf,meancorr_uncf),[],{'bin','confidence'});disp(tbl);
fprintf('eta_p_squared = %.5f\n\n',tbl.SumSq(5)/sum(tbl.SumSq(5:6)))
%
X = [meanconf_cor(:,4) mean(meanconf_cor(:,[3 5]),2) mean(meanconf_cor(:,[2 6]),2) mean(meanconf_cor(:,[1 7]),2)];
Y = [meanconf_err(:,4) mean(meanconf_err(:,[3 5]),2) mean(meanconf_err(:,[2 6]),2) mean(meanconf_err(:,[1 7]),2)];
fprintf('ANOVA : confidence explained by stimulus strength on CORRECT trials\n');
tbl = simple_mixed_anova(X,[],{'strength'}); disp(tbl);
fprintf('eta_squared = %.5f\n\n',tbl.SumSq(3)/sum(tbl.SumSq(3:4)))
%
fprintf('ANOVA : confidence explained by stimulus strength on ERROR trials\n');
tbl = simple_mixed_anova(Y,[],{'strength'}); disp(tbl);
fprintf('eta_squared = %.5f\n\n',tbl.SumSq(3)/sum(tbl.SumSq(3:4)))

% plot fraction correct and relative confidence depending on luminance bins
pbar = 2.8; HGT = 2.5; xvec = -3:+3;
% plot fraction correct
yavg_g = mean(meancorr_conf,1);
yerr_g = std(meancorr_conf,[],1)/sqrt(ndata);
yavg_r = mean(meancorr_uncf,1);
yerr_r = std(meancorr_uncf,[],1)/sqrt(ndata);

figure('Color','white','Name','Figure_5c_luminance_accuracy');

hold on
ylim([0.5 1]);
patch([xvec,fliplr(xvec)],[yavg_g+yerr_g,fliplr(yavg_g-yerr_g)], ...
    0.5*([0 1 0]+1),'EdgeColor','none');
patch([xvec,fliplr(xvec)],[yavg_r+yerr_r,fliplr(yavg_r-yerr_r)], ...
    0.5*([1 0 0]+1),'EdgeColor','none');
for i = -3.5:3.5
    plot([i,i],ylim,'k:');
end
plot(xvec,yavg_g,'Color','g','LineWidth',1.5);
plot(xvec,yavg_r,'Color','r','LineWidth',1.5);
set(gca,'XTick',-2.5:2.5,'YTick',0.5:0.1:1);
xlabel('marble luminance'); ylabel('fraction correct (%%)');
helper_ILL(gcf,pbar,HGT);

% plot relative confidence
yavg_g = mean(meanconf_cor,1);
yerr_g = std(meanconf_cor,[],1)/sqrt(ndata);
yavg_r = mean(meanconf_err,1);
yerr_r = std(meanconf_err,[],1)/sqrt(ndata);

figure('Color','white','Name','Figure_5c_luminance_confidence');
hold on
ylim([-.5 .5]);
patch([xvec,fliplr(xvec)],[yavg_g+yerr_g,fliplr(yavg_g-yerr_g)], ...
    0.5*([.5 .5 .5]+1),'EdgeColor','none');
patch([xvec,fliplr(xvec)],[yavg_r+yerr_r,fliplr(yavg_r-yerr_r)], ...
    0.5*([.5 .5 .5]+1),'EdgeColor','none');
for i = -3.5:3.5
    plot([i,i],ylim,'k:');
end
plot(xvec,yavg_g,'Color',[.5 .5 .5],'LineWidth',1.5);
plot(xvec,yavg_r,'Color',[.5 .5 .5],'LineWidth',1.5,'LineStyle',':');
set(gca,'YTick',-1:0.5:1);
xlabel('marble luminance'); ylabel('relative confidence');
helper_ILL(gcf,pbar,HGT);

%% Characterizing the suboptimality of statistical inference (replication from experiment 1)
clearvars;
expe2_analysis_suboptimality

%% Response stabilization mechanisms (predictions validation)
clearvars;
expe2_analysis_stabilization