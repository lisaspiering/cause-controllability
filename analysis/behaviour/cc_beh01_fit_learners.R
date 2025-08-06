# Script to fit the 3 different Bayesian learning models
# 1) Active learner
# 2) Passive learner
# 3) Ignorant loearner

rm(list=ls()) 
library(tidyverse);library(dplyr);library(ggplot2);library(ggpubr);library(rstan);library(brms);library(shinystan);library(rmcorr);library(zoo);library(foreach);library(doParallel);library(sjPlot);library(cmdstanr)
registerDoParallel(cores=14)# for parallel computing 

# fitting variables
learners2fitall = c('BL','IL','PL') # BL = Active learner, IL = Ignorant learner, PL = Passive learner
sample2fitall   = c('mri','online')
phase2fitall    = c(1,3)        # 1 or 3, or both; 1 = Self-Other phase, 3 = Control-Other phase
sd2fit          = 30 # initial SD = uncertainty
max.tries       = 3
last_pre_drift_trial = 11 # trial number that is the last one pre-drift

### Source  functions ###
sapply(list.files(path='/Users/lisaspiering/Documents/GitHub/cause-controllability/analysis/_functions/', pattern='.R', full.names = T),source) 

# load MRI, online & BL data
out  = setup_cc()
path = out$path; subs2filter = out$subs2filter; df2plot = out$df2plot; beta2plot = out$beta2plot; my.priors = out$my.priors
setwd(path$script) # set working directory to script folder so that stan finds stan file

# load data
# 1) load MRI behaviour
load(file = file.path(path$data,'preprocessed_behaviour_mri.RData'),verbose=T)
dataB.mri = dataB.mri %>% filter(ID %in% subs2filter$mri.incl) 

# 2) load online data
load(file.path(path$data,'preprocessed_behaviour_online.RData'),verbose=T)
dataB.online = dataB.online %>% filter(!(ID %in% subs2filter$online.excl))

# for our regressions, I will look at both online and MRI participants at once -> for this, combine both datasets in one and therefore create a new ID column
dataB.all = data.frame(); i = 1
for (s in subs2filter$mri.incl) {
  dataB.all = rbind(dataB.all,
                    dataB.mri %>% filter(ID==s,phase!=2) %>% select(ID,game_number,phase,phase_trial,player_mu_drift,feedback,feedback_mu,self_mu, player_mu,control_level,active_trial,rate_self,rate_other,rate_ctrl) %>%
                      mutate(ID2 = i,sample='mri'))
  i = i+1
}
for (s in unique(dataB.online$ID)) {
  dataB.all = rbind(dataB.all,
                    dataB.online %>% filter(ID==s,phase!=2) %>% select(ID,game_number,phase,phase_trial,player_mu_drift,feedback,feedback_mu,self_mu, player_mu,control_level,active_trial,rate_self,rate_other,rate_ctrl) %>%
                      mutate(ID2 = i,sample='online'))
  i = i+1
}

dataB.all.BL = dataB.all
# set up ignorant learner: here we can just set active_trial = 0 
dataB.all.IL = dataB.all %>% mutate(active_trial=0)
# set up passive learner : this one doesn't learn from active trial -> I remove those trials
dataB.all.PL = dataB.all %>% 
  # for my fitting, I need people's priors on first trial & the last trial before drift -> if that one was an active trial, I just replace it by a random feedback around the feedback_mu
  rowwise() %>% mutate(feedback = ifelse(phase_trial%in%c(1,11) & active_trial==1, rnorm(1,mean=feedback_mu,sd=5), feedback),
                       active_trial = ifelse(phase_trial%in%c(1,11) & active_trial==1, 0, active_trial)) %>% # on first phase_trial, set active_trial to 0 so that we can extract the priors here when fitting for the PL
  filter(active_trial==0) 

# start fitting
for (learner2fit in learners2fitall) {
  for (phase2fit in phase2fitall){
    ## take respective sample's data & set path where to save to
    if(learner2fit=='BL') {data2fit=dataB.all.BL} else if(learner2fit=='IL') {data2fit=dataB.all.IL}  else if(learner2fit=='PL') {data2fit=dataB.all.PL} else { warning('learner2fit unknown')}
    path$save2 = file.path(path$data,learner2fit)
    
    # depending on phase, remove columns we won't need
    if (phase2fit==1) {
      data2fit     = data2fit %>% select(ID,ID2,game_number,phase,player_mu_drift,phase_trial,active_trial,rate_self,rate_other,control_level,feedback,sample) %>% filter(phase==phase2fit)
      stanfile     = 'BayesianModel_phase1.stan'
      pars2extract = c("RateSelfMean", "RateOtherMean","RateOtherMeanDrift")
      priors       = list(names  = c('PriorSelfMean','PriorSelfSD','PriorOtherMean','PriorOtherSD','PriorOtherMeanDrift'), 
                          values = c('rate_self',     sd2fit,      'rate_other',     sd2fit,       'RateOtherMean')) 
    } else if (phase2fit==3){
      data2fit     = data2fit %>% select(ID,ID2,game_number,phase,player_mu_drift,phase_trial,active_trial,rate_self,rate_other,rate_ctrl,feedback,self_mu,sample) %>% mutate(self_sd=5)  %>% filter(phase==phase2fit)
      stanfile     = 'BayesianModel_phase3.stan'
      pars2extract = c("RateOtherMean","RateControlMean","RateOtherMeanDrift")
      priors       = list(names  = c('PriorOtherMean','PriorOtherSD','PriorControlMean','PriorControlSD','PriorOtherMeanDrift'), 
                          values = c('rate_other',     sd2fit,           'rate_ctrl',     sd2fit/100,     'RateOtherMean')) 
    }
    
    # precompile the model
    mod = cmdstan_model(stan_file = stanfile)
    
    # loop over participants
    foreach(s=sort(unique(data2fit$ID2))) %dopar% {
      # set up some empty dataframes to store results in
      bay_means.df = data.frame(); bay_conf_up.df = data.frame(); bay_conf_lo.df = data.frame(); bay_fit.df = data.frame()
      idata      = data2fit %>% filter(ID2==s) 
      sample2fit = idata %>% pull(sample) %>% unique()
      s_orig     = idata %>% pull(ID) %>% unique()# get original participant ID for saving this
      # loop over games
      for (g in sort(unique(idata$game_number))){
        # loop over trials
        for (itrial in sort(unique(idata$phase_trial))){
          # extract schedule data up to this trial and add prior beliefs to it
          df          = idata %>% filter(game_number==g,phase_trial<=itrial) %>% select(-sample,-ID,-ID2)
          if(nrow(df)>0){
            df          = lapply(as.list(df),as.array) # otherwise Stan treats it as a scalar, rather than as array
            df$nrow     = length(df$feedback)
            df[[priors$names[1]]] = idata %>% filter(game_number==g,phase_trial==1) %>% pull(priors$values[1]) # take participants' prior from first trial's rating
            df[[priors$names[2]]] = as.numeric(priors$values[2]) # standard deviation
            df[[priors$names[3]]] = idata %>% filter(game_number==g,phase_trial==1) %>% pull(priors$values[3]) # take participants' prior from first trial's rating
            df[[priors$names[4]]] = as.numeric(priors$values[4]) # standard deviation
            # set prior of other player post-drift as the last estimated value from the learner pre-drift (is saved in bay_means.df)
            df[[priors$names[5]]] = ifelse(itrial>last_pre_drift_trial, 
                                           bay_means.df%>%filter(game_number==g,phase_trial==last_pre_drift_trial)%>%pull(priors$values[5]), # as prior for post-drift, take learner's last estimate 
                                           50) # 50 for all pre-drift trials (but doesn't matter, we don't use 50 really)
            
            # initialise fitting parameters to keep track of fitting to criterion
            count_divergent = 1 # pseudo values so that while loop below works
            highest_rhat    = 5 # pseudo values so that while loop below works
            my.iter         = 10000
            try.count       = 1
            
            # run stan
            while (count_divergent > 0 || highest_rhat >= 1.1){
              suppressWarnings(rm(simDatT)) # rm output if it exists, and suppress warning if it doesn't exist yet
              simDatT <- mod$sample(data=df,chains=4,iter_sampling=my.iter, adapt_delta = 0.99,refresh=0) 
              
              # get fitting things
              count_divergent  = simDatT$diagnostic_summary(diagnostics = 'divergences') %>% unlist() %>% sum()
              highest_rhat    = rhat(simDatT)[names(rhat(simDatT))!=c('lp__')] %>% max() # get highest rhat value (don't check rhat for lp__)
              my.iter         = my.iter*1.5
              try.count       = try.count+1
              if (try.count>max.tries){
                break
              }
            }
            
            # summarise the results
            bay_mean    = simDatT$summary(variables = pars2extract) %>% select(variable,mean) %>% pivot_wider(names_from = variable,values_from = mean) %>% as.data.frame()
            bay_conf_up = simDatT$summary(variables = pars2extract) %>% select(variable,q95)  %>% pivot_wider(names_from = variable,values_from = q95)  %>% as.data.frame()
            bay_conf_lo = simDatT$summary(variables = pars2extract) %>% select(variable,q5)   %>% pivot_wider(names_from = variable,values_from = q5)   %>% as.data.frame()
            bay_fit     = simDatT$summary(variables = pars2extract) %>% select(variable,rhat) %>% pivot_wider(names_from = variable,values_from = rhat) %>% as.data.frame()
            
            # extract the relevant values
            tDat1       = cbind(bay_mean,phase_trial=itrial,game_number=g,ID2=s,ID=s_orig,sample=sample2fit)
            tDat2       = cbind(bay_conf_up,phase_trial=itrial,game_number=g,ID2=s,ID=s_orig,sample=sample2fit)
            tDat3       = cbind(bay_conf_lo,phase_trial=itrial,game_number=g,ID2=s,ID=s_orig,sample=sample2fit)
            tDat4       = cbind(bay_fit,phase_trial=itrial,game_number=g,ID2=s,ID=s_orig,sample=sample2fit) %>% mutate(iter = my.iter, div_transitions = count_divergent, fitted.ok.flag = ifelse(try.count>max.tries, 0, 1))
            
            # bind the relevant values into overall dataframes
            bay_means.df     = rbind(bay_means.df,tDat1)
            bay_conf_up.df   = rbind(bay_conf_up.df,tDat2)
            bay_conf_lo.df   = rbind(bay_conf_lo.df,tDat3)
            bay_fit.df       = rbind(bay_fit.df,tDat4) 
            
            # clean up
            rm(simDatT,df,tDat1,tDat2,tDat3,bay_mean,bay_conf_up,bay_conf_lo,bay_fit)
          } else {warning('Empty data!')}
        }
      }
      # save simulation results into one variable and write to disk
      sim_results = list('Means' = bay_means.df, 'Upper confidence' = bay_conf_up.df, 'Lower confidence' = bay_conf_lo.df, 'Fit' = bay_fit.df, 'idata' = idata, 'sd2fit' = sd2fit, 'learner2fit'=learner2fit)
      # save
      fileName   = paste0(learner2fit,'_p',phase2fit,'_s',str_pad(s_orig,2,pad='0'),'.RData')
      save(sim_results, file = file.path(path$save2,sample2fit,fileName))
    }
  }
}

