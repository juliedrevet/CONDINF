%% RUNNER_ANALYSIS
% needed on the path for running analysis:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)
%       * simple_mixed_anova.m (Calpette, L. (2022) https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)
clearvars
close all
clc
addpath '../helper_plot';
addpath '../matlab_functions';
datapath = '../data_table';

%% look for outliers depending on overall performance across both experiments
% load behavioural data
nstd = 3.5; 
pcor = [];
subjlist = [];
for iexpe = 1:2
    % load data table
    switch iexpe
        case 1
            load(fullfile(datapath,'expe1_data.mat'),'expe1_data');
            dat = expe1_data; clear('expe1_data');
            participants = unique(dat.participant);
            ndata = numel(participants);
            subjlist = cat(1,subjlist,[ones(ndata,1) participants]);
        case 2
            % load data table
            load(fullfile(datapath,'expe2_data.mat'),'expe2_data');
            dat = expe2_data; clear('expe2_data');
            participants = unique(dat.participant);
            ndata = numel(participants);
            subjlist = cat(1,subjlist,[2*ones(ndata,1)  participants]);
    end
    
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata);
        s = dat.stim(ifilt);    % stimulus
        r = dat.resp(ifilt);    % response
        % compute performance
        pcor = cat(1,pcor,nanmean(r == s));
    end
end
% compute leave-one out group-level mean performance
ndata = size(pcor,1);
pcor_loo_avg = nan(ndata,1);
pcor_loo_std = nan(ndata,1);
for idata = 1:ndata
    pcor_loo_avg(idata) = mean(pcor(setdiff(1:ndata,idata)));
    pcor_loo_std(idata) = std(pcor(setdiff(1:ndata,idata)));
end
I = find((pcor<(pcor_loo_avg+pcor_loo_std*(-nstd))) | (pcor > (pcor_loo_avg+pcor_loo_std*(+nstd))));
if ~isempty(I)
    fprintf('remove participant S%02d from analysis (Experiment %d)\n',subjlist(I,2),subjlist(I,1))
end
%% run all analyses 
runner_analysis_expe1;

runner_analysis_expe2;

runner_analysis_expe_merged;

runner_analysis_supplementary;

