# controllability and cause in human collaboration

Spiering et al., 2024.
bioRxiv preprint available at: https://www.biorxiv.org/content/10.1101/2024.10.01.615833v1

This repository contains data, analysis script and plotting scripts for the above project. Please refer to above preprint for details on the project. Each script contains comments on how to use the script and which figures or analyses it produces. For any questions, please refer to above preprint or contact the authors.

## system requirements

The code was developed on macOS Sonoma 14.7.6, using R version 4.2.1 (2022-06-23), FSL (v. 6.00) and MATLAB (R2021a). A list of required R packages can be found in the `R-packages.txt` file

## installation guide

Please install these programs:

- R: https://www.r-project.org/
- FSL: https://fsl.fmrib.ox.ac.uk/fsl/docs/#/install/index
- MATLAB: https://uk.mathworks.com/help/install/ug/install-products-with-internet-connection.html

Then install the required R packages (see R-packages.txt) using the `install.packages()` function in R. This can take a few minutes.

## instruction for use

**analysis:**

The analysis folder contains the scripts necessary to analyse the behavioura and neural data provided, and reproduce the plots of the study.

- `_functions` : contains custom R functions used for the analysis
  - `setup_cc.R`: setup the analysis environment, load the required packages and data. change the paths here to run the analyses and plotting scripts
  - `fit_sub_behaviour.R`: fit subject-wise regressions using brms
  - `load_data4regs.R`: load the preprocessed behavioural data to run regression analyses and for plotting
  - `load_learners.R`: load the simulated data from the three Bayesian optimal observer models and aggregate with the behavioural data
  - `load_tc2plot.R`: load the time course data for plotting the fMRI time courses
  - `transform_data_l1.R`: transform the behavioural data in a long format for the fMRI time course analyses
- `behaviour`: contains analysis and plotting script for the behavioural results
  - `BayesianModel_phase1.stan` and `BayesianModel_phase3.stan`: the stan files that specify the Bayesian learner models (phase1 = Self-Other phase, phase3 = Control-Other phase).
  - `cc_beh01_fit_learners.Rmd`: fits the Bayesian learner models to the behavioural data
  - `cc_beh02_overview_behaviour.Rmd`: plotting script for various overview plots of the behavioural data
  - `cc_beh03_credit_assignment_social_PEs.Rmd`: regression analysis and plotting script for how participants assign social PEs to themselves and others according to their perceived control (Figure 2e)
  - `cc_beh04_AD.Rmd`: regression analysis and plotting script for how participants used AD (e.g. Figure 3)
  - `cc_beh05_regression_intercept_anovas.Rmd`: analysis script to test the intercept effects of the resulting beta estimates from subject-wise regressions (e.g. Figure 2e)
  - `cc_beh06_AD_guided_control_learning.Rmd`: regression analysis and plotting script for how participants used AD to guide their control learning (e.g. Figure 4)
- `fmri`: contains analysis and plotting script for the fMRI results 
  - subfolder `functions` contains the custom MATLAB functions supporting the script `cc_roi02_run_timecourse.mlx`
  - `cc_roi01_prep_behavioural_data.Rmd`: prepares the behavioural data for the fMRI time course analyses
  - `cc_roi02_run_timecourse.mlx`: runs the fMRI time course analyses and performs LOO-procedure
  - `cc_roi03_plot.Rmd`: plots the fMRI time courses (e.g. Figures 6 and 7)

**data:**

The behavioural data can be found in the folder `behavioural-data` and the fMRI data can be found in the folder `mri-data`.

- `behavioural-data`: 
  - The behavioural data for the respective behavioural analyses are found in the files `preprocessed_behaviour_mri.RData` (mri sample), `preprocessed_behaviour_online.RData` (online sample), and `preprocessed_data_4regs.RData` (both samples combined, for the regressions)
  - The regression fits are found in the subfolder `regression-fits`. More information on which files are necessary for which analysis or plotting script can be found in the behavioural analysis scripts folder (`analysis\behaviour`)
- `mri-data`:
  - The cluster-corrected whole-brain maps shown in Fig 5, showing the effects of uncertainty and AD at the time of action, can be found in the folder `whole-brain-maps` in `mri-data`. The cluster-corrected maps are denoted as `*_thresh_act.nii.gz` and `_thresh_deact.nii.gz`, and the uncorrected maps are denoted as `*_unthresholded.nii.gz`. These can be viewed using for example FSLeyes.
  - The ROI masks in MNI space used for the fMRI time-course analyses can be found in the folder `roi-masks` in `mri-data`.
  - For the fMRI time-course analysis, the time-course data for each ROI can be found in the folder `roi-timecourses` in `mri-data`. The files are named according to to the script they were generated with, e.g. `roi01_behdata_glm001.mat` contains the behavioural regressor data for the ROI-GLMs defined in the script `cc_roi01_prep_behavioural_data.Rmd`.

**demo** : The code can be used together with the data provided to reproduce various figures and analyses for the above project (see above).

Copyright, 2025, Lisa Spiering. All rights reserved
