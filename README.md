# CONDITIONAL BELIEF UPDATING

This is the overall repository that replicates all analyses from
Drevet, J., Drugowitsch, J. & Wyart, V. Efficient stabilization of imprecise statistical inference through conditional belief updating. _Nat Hum Behav_ (2022). https://doi.org/10.1038/s41562-022-01445-0

Participants' behavioral data can be found in the repository data_table or on figshare:
https://figshare.com/projects/Efficient_stabilization_of_imprecise_statistical_inference_through_conditional_belief_updating/140170

To fit the models, you need to add following code to your path:
- VBMC v1.0 (https://github.com/lacerbi/vbmc/)
- BADS v1.0 (https://github.com/lacerbi/bads/)
      
To run the analyses, you need to add following functions to your path:
- spm_BMS.m in SPM12 (Wellcome Trust Center for Human Neuroimaging; http://www.fil.ion.ucl.ac.uk/spm)
- simple_mixed_anova.m (Calpette, L. (2022) https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)
    
### This repository contains ###
- data_table: contains participants' behavioral data as .mat
- results: contains all model-fitting results as .mat
- runner_analysis: contains all Matlab runner scripts to replicate the analyses and figures of the paper.
- runner_fit: contains all Matlab runner scripts to replicate the model-fitting. All fits are already available in the results folder.
- helper_plot: contains scripts to config plots. Automatically added to path in the runner_analysis scripts.
- matlab_functions: contains all scripts used by runner scripts. Automatically added to path in the runner_analysis or runner_fit scripts.

### What to do to replicate the analyses and figures from the paper ###
- make sure you have downloaded following folders: data_table, runner_analysis, results, matlab_functions, helper_plot
- go to the runner_analysis folder
- to run ALL analyses from the paper (including supplementary): RUNNER_ANALYSIS.m
- you can also run individual scripts:

      Figures 1-4: runner_analysis_expe1.m
      Figures 5-6: runner_analysis_expe2.m 
      Figures 7-8: runner_analysis_expe_merged.m
      all Supplementary Figures: runner_supplementary.m

     
### What to do to replicate the model-fitting results (LONG!!) ###
All model-fitting results are already stored in the results folder. To run them again (note that this takes a lot of time):
- make sure you have downloaded following folders: data_table, matlab_functions, runner_fit
- got to the runner_fit folder
- For experiment 1: Run the script runner_fit_expe1.m 
- For experiment 2: Run the script runner_fit_expe2.m 
