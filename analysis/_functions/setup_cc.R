setup_cc <- function() {
  # this function sets some general variables relevant for the data analysis
  # we also set some settings for the plotting (e.g. the general ggplot theme) and source all subfunctions
  # inputs:
  #  - paths of the data
  # outputs:
  #  - subjects2filter: list of participant IDs to in- or exclude in this study
  #  - path: list of paths that are relevant for this study
  #  - my.priors: variable that sets for brms regressions the priors for the beta estimates
  #  - df2plot + beta2plot: 2 plotting variables that help us, for regressions, to plot either the normal beta estimates or the log beta estimates
  
  # set seed
  set.seed(123)
  
  # define paths
  path        = list(base = file.path('/Users',Sys.info()['user'],'Documents'))
  path$local  = file.path('/Volumes/G_DRIVE_Lisa','mood-fmri') 
  path$gitlab = file.path(path$base,'GitHub','cause-controllability')
  path$onedrive = file.path('/Users',Sys.info()['user'],'OneDrive - Nexus365/Oxford/_03_PhD MRI study/')
  path$save   = file.path(path$local,'_behaviour','brms_healthy')
  path$data   = file.path(path$gitlab,'behavioural-data')
  path$script = file.path(path$gitlab,'analysis','behaviour')
  path$figs   = file.path(path$onedrive,'_manuscript/figures')
  path$tables = file.path(path$onedrive,'_manuscript/tables')
  
  # make a list of participants to exclude and include (mri and online sample)
  subs2filter = list(
    # MRI participants that we want to include
    mri.incl    = c(1,  3,  4,  5,  6,  8,  9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34),
    # ONLINE participants to EXclude (who don't have enough AD trials in their behaviour, reported never having done AD, or who rated themselves lower after/before AD)
    online.excl = c(1,  5,  6, 13, 15, 19, 20, 21, 23, 24, 25, 26, 31, 33, 36, 37, 39, 40, 47, 52, 53, 55, 59, 60, 61, 62, 63, 64, 66, 67))
  
  # brms prior settings
  my.priors = set_prior("normal(0, 5)", class = "b")
  
  # set theme for all ggplots
  theme_set(theme_pubr(base_size=14))
  
  # plotting variables
  log2plot = 0 # for regressions, do we want to plot the beta estimates or the log transformed beta estimates? -> we do non-transformed beta estimates
  df2plot = ifelse(log2plot,'df.log','df');  beta2plot = ifelse(log2plot,'log_beta_estimate','beta_estimate') 
  cols = list('Control'  = '#CC79A7', 'Self' = '#009E73', 'Other' = '#0072B2', 'AD' = '#E69F00',
              'ControlL' = '#F7D0E0', 'SelfL'= '#B8E8D0', 'OtherL' = '#AFD0F1','ADL' = '#FFDC9A')
  
  # give these variables back
  out = list(path = path, subs2filter = subs2filter, df2plot = df2plot, beta2plot = beta2plot, my.priors = my.priors,cols=cols)
  return(out)
}