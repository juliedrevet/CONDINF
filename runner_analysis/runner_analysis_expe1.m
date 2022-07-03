%% RUNNER ANALYSIS EXPERIMENT 1 
% needed on the path for running analysis:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)

addpath '../helper_plot';
addpath '../matlab_functions';

fprintf('<strong>***\tEXPERIMENT 1\t***</strong>\n\n')

%% Psychophysical titration procedure
fprintf('\n<strong>Psychophysical titration procedure</strong>\n\n')

clearvars;
datapath = '../data_table/';
load(fullfile(datapath,'expe1_titration.mat'));
dat = expe1_titration;
clearvars expe1_titration;
participants = setdiff(unique(dat.participant),33); % participant S33 removed from analysis
ndata = numel(participants);
nblck = numel(unique(dat.block));

% get luminance thresholds per participant
xmin = zeros(ndata,nblck);
xmax = zeros(ndata,nblck);

for idata = 1:ndata
    ifilt = dat.participant == participants(idata) & dat.trial == 1;
    xmin(idata,:) = dat.xmin(ifilt);
    xmax(idata,:) = dat.xmax(ifilt);
end
% averaging over blocks before taking median over participants
xmin_median = median(mean(xmin));
xmax_median = median(mean(xmax));
dark_grey  = mean(xmin,2);
light_grey = mean(xmax,2);

fprintf('Dark grey color ratio  (mean +/- s.e.m.): %4.4f +/- %4.4f\n',mean(dark_grey),std(dark_grey)/sqrt(ndata))
fprintf('Light grey color ratio (mean +/- s.e.m.): %4.4f +/- %4.4f\n',mean(light_grey),std(light_grey)/sqrt(ndata))
fprintf('Dark grey color ratio (median)  : %4.4f\n',mean(xmin_median))
fprintf('Light grey color ratio (median) : %4.4f\n\n',mean(xmax_median))

% psychometric curves
fprintf('<strong>Computing regression for psychometric curves\n</strong>')
fun = @(pp,xx) pp(:,3)/2+normcdf(pp(:,1)+pp(:,2)*xx)*(1-pp(:,3));
xfun = 0:.01:1;
yfun = zeros(nblck,length(xfun),ndata);
bfun = zeros(nblck,3,ndata);

c = lines;
w = 0.5;

for idata = 1:ndata
    fprintf('\tParticipant %02d/%02d.',idata,ndata)
    for iblck = 1:nblck
        fprintf('.')
        ifilt = dat.participant == participants(idata) & dat.block <= iblck & dat.session == unique(dat.session(dat.block == iblck));
        xdata = dat.prop(ifilt);
        ydata = dat.resp(ifilt) == 1;
        
        wdata = zeros(sum(ifilt),1);
        for jblck = 1:iblck
            wdata(dat.block(ifilt) == jblck) = w.^(mod(iblck-jblck,8));
        end

        % perform logistic regression on categorized ratio
        b = logreg(fun,xdata,ydata,wdata);
        bfun(iblck,:,idata) = b;
        yfun(iblck,:,idata) = fun(b,xfun);

    end
    fprintf('done.\n')
end

dark_cues  = dat.prop(dat.staircase<0);
light_cues = dat.prop(dat.staircase>0);

% plot
figure('Color','white','Name','Figure_1c_titration_procedure');
shft = .4; pbar = .9; HGT = 5;
rgb = [185 208 84]/255;
hold on
xlim([.35 .65])
h1 = histogram(dark_cues,20,'FaceColor',[66,66,66]/255,'EdgeColor','none','Normalization','probability','FaceAlpha',1);
h2 = histogram(light_cues,20,'FaceColor',[190,190,190]/255,'EdgeColor','none','Normalization','probability','FaceAlpha',1);
fill([xmin_median xmax_median xmax_median xmin_median],[0 0 1 1]+shft,[.8 .8 .8],'EdgeColor','none');
for idata = 1:ndata
    p = plot(xfun,mean(yfun(:,:,idata),1)+shft,'Color',(rgb+1)/2,'LineWidth',.5);
end
% shaded area = median threshold level across all participants.
plot(xfun,mean(mean(yfun,3),1)+shft,'Color',rgb,'LineWidth',1.5)
line([0 xmin_median],[.2+shft .2+shft],'LineStyle',':','Color','k')
line([0 xmax_median],[.8+shft .8+shft],'LineStyle',':','Color','k')
line([xmin_median xmin_median],[shft .2+shft],'LineStyle',':','Color','k')
line([xmax_median xmax_median],[shft .8+shft],'LineStyle',':','Color','k')
xlim([.35 .65])
ylabel('fraction perceived as light'); xlabel('marble luminance')
set(gca,'XTick',0:.1:1,'YTick',[.2 .8]+shft,'YTickLabel',{'20%', '80%',});
helper_ILL(gcf,pbar,HGT);

%% Characterizing the suboptimality of statistical inference
clearvars;
expe1_analysis_suboptimality

%% Response stabilization mechanisms
clearvars;
expe1_analysis_stabilization