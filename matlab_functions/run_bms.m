function [pavg,pstd,pexc,alpha] = run_bms(lme)
%  RUN_BMS  Run Bayesian model selection using Dirichlet parameterization
% needed on the path:
%       * spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)

% check presence of spm_BMS function
if ~exist('spm_BMS','file')
    error('Missing spm_BMS function in path!');
end

% run Bayesian model selection using Dirichlet parameterization
[alpha,~,pexc] = spm_BMS(lme);

% get means and variances for the Dirichlet parameterization
[pavg,pvar] = drchstat(alpha);
pstd = sqrt(pvar);

end