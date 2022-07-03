# CONDITIONAL BELIEF UPDATING

This is the overall repository that replicates all analyses from Drevet J, Drugowitsch J, Wyart V (2022) _Efficient stabilization of imprecise statistical inference through conditional belief updating_ https://www.biorxiv.org/content/10.1101/2022.06.08.495322v1

Participants' Data can be found on figshare:
https://figshare.com/projects/Efficient_stabilization_of_imprecise_statistical_inference_through_conditional_belief_updating/140170

To fit the models, you need to add following code to your path:
- VBMC v1.0 (https://github.com/lacerbi/vbmc/)
- BADS v1.0 (https://github.com/lacerbi/bads/)
      
To run the analyses, you need to add following functions to your path:
- spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)
- simple_mixed_anova.m (Calpette, L. (2022) https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)
    
### This repository contains ###
- results: contains all model-fitting results as .mat
- runner_analysis: contains all Matlab runner scripts to replicate the analyses and figures of the paper.
- runner_fit: contains all Matlab runner scripts to replicate the model-fitting. All fits are already available in the results folder.
- helper_plot: contains scripts to config plots. Automatically added to path in the runner_analysis scripts.
- matlab_functions: contains all scripts used by runner scripts. Automatically added to path in the runner_analysis or runner_fit scripts.

### What to do to replicate the analyses and figures from the paper ###
- make sure you have downloaded following folders: runner_analysis, results, helper_plot
- go to the runner_analysis folder
- to run ALL analyses from the paper (including supplementary): RUNNER_ANALYSIS.m
- you can also run individual scripts for experiment 1 (runner_analysis_expe1.m), for experiment 2 (runner_analysis_expe2.m), for analyses based on both experiments merged (runner_analysis_expe_merged.m), or for supplementary information (runner_supplementary.m)
     
### What to do to replicate the model-fitting results (LONG!!) ###
All model-fitting results are already stored in the results folder. To run them again (note that this take a lot of time):
- download Participants' Data from figshare and save them in a folder called data_table
- For experiment 1: Run the script runner_fit_expe1.m
- For experiment 2: Run the script runner_fit_expe2.m
