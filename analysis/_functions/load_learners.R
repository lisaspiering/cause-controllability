load_learners = function(dataB.mri,dataB.online,path) {
  # this function is for loading the simulated data from the Bayesian learning models
  
  # loop through the subfolders in the saved BL results  / mri and online
  for (learner2load in c('BL','IL','PL')) { # PL should be last so that the left_join below works best (bc PL doesn't have all data rows as IL & BL)
    
    for (ifolder in c('mri','online')) {
      
      suppressWarnings(rm(dataP))
      if (ifolder=='mri') {assign('dataP',dataB.mri)} else {assign('dataP',dataB.online)} # copy participants' data (mri or online) 
      
      # get all the fit BL files 
      allBL  = list.files(file.path(path$data,learner2load,ifolder),full.names = T,pattern = '.RData')
      
      # predefine an empty dataframe where we want to store all the BL results in (BL_ prefix will be replaced later with correct learner2load prefix)
      colnam = c("ID2","sample","ID","game_number","phase_trial","active_trial","rate_self","rate_other","rate_ctrl","feedback","self_mu","self_sd","phase",
                 "BL_self_up","BL_self_mu","BL_self_lo","BL_other_up","BL_other_mu","BL_other_lo","BL_ctrl_up","BL_ctrl_mu","BL_ctrl_lo",
                 'BL_self_rhat','BL_other_rhat','BL_ctrl_rhat','BL_niter','BL_ndiv','BL_fittedok')
      
      dataBL = data.frame(matrix(nrow=0,ncol=length(colnam))); colnames(dataBL) <- colnam
      
      # now loop through all the BL fits and put it together into dataBL
      for (ifile in allBL) {
        
        load(ifile) # load the file
        iphase = str_split(ifile,'_')[[1]][2] # get the phase; 'p1' or 'p3'
        # depending on the phase, extract the relevant BL estimations (self + other for phase 1, other + ctrl for phase 3)
        if (iphase=='p1') {
          sim_results$idata              = sim_results$idata              %>% mutate(phase=1,rate_ctrl=NA,self_mu=NA,self_sd=NA,BL_ctrl_mu=NA,BL_ctrl_up=NA,BL_ctrl_lo=NA,BL_ctrl_rhat=NA) # need to add some empty columns so that we can put p1 and p3 together into one dataframe
          sim_results$Means              = sim_results$Means              %>% rename(BL_self_mu = RateSelfMean, BL_otherpre_mu = RateOtherMean, BL_otherpost_mu = RateOtherMeanDrift) 
          sim_results$`Upper confidence` = sim_results$`Upper confidence` %>% rename(BL_self_up = RateSelfMean, BL_otherpre_up = RateOtherMean, BL_otherpost_up = RateOtherMeanDrift)
          sim_results$`Lower confidence` = sim_results$`Lower confidence` %>% rename(BL_self_lo = RateSelfMean, BL_otherpre_lo = RateOtherMean, BL_otherpost_lo = RateOtherMeanDrift)
          # store fitting diagnostics from Bayesian Learner
          sim_results$Fit                = sim_results$Fit %>% rename(BL_self_rhat = RateSelfMean, BL_otherpre_rhat = RateOtherMean,BL_otherpost_rhat = RateOtherMeanDrift, BL_niter = iter, BL_ndiv = div_transitions, BL_fittedok = fitted.ok.flag)
          
        } else {
          sim_results$idata              = sim_results$idata              %>% mutate(phase=3,control_level=NA,BL_self_mu=NA,BL_self_up=NA,BL_self_lo=NA,BL_self_rhat=NA)
          sim_results$Means              = sim_results$Means              %>% rename(BL_ctrl_mu = RateControlMean, BL_otherpre_mu = RateOtherMean, BL_otherpost_mu = RateOtherMeanDrift)
          sim_results$`Upper confidence` = sim_results$`Upper confidence` %>% rename(BL_ctrl_up = RateControlMean, BL_otherpre_up = RateOtherMean, BL_otherpost_up = RateOtherMeanDrift)
          sim_results$`Lower confidence` = sim_results$`Lower confidence` %>% rename(BL_ctrl_lo = RateControlMean, BL_otherpre_lo = RateOtherMean, BL_otherpost_lo = RateOtherMeanDrift)
          sim_results$Fit                = sim_results$Fit                %>% rename(BL_otherpre_rhat = RateOtherMean, BL_otherpost_rhat = RateOtherMeanDrift, BL_ctrl_rhat = RateControlMean, BL_niter = iter, BL_ndiv = div_transitions, BL_fittedok = fitted.ok.flag)
        }
        
        # put mean, upper confidence, lower confidence results from list format together in one data frame
        isim = left_join(sim_results$idata,sim_results$`Upper confidence`,by=c('ID2','sample','ID','game_number','phase_trial')) %>% 
          left_join(sim_results$Means,by=c('ID2','sample','ID','game_number','phase_trial')) %>% 
          left_join(sim_results$`Lower confidence`,by=c('ID2','sample','ID','game_number','phase_trial')) %>%
          left_join(sim_results$Fit, by=c('ID2','sample','ID','game_number','phase_trial')) # add fitting diagnostics
        
        # put this BL in big dataframe with everybody else
        dataBL = rbind(dataBL,isim)
      }
      
      # add BL to rest of the mri participants' data - for this we need to keep in mind that we need to shift BL ratings by one trial (because BL estimation of observing feedback of trial 1 are saved in trial 1 / not in trial 2 like for participants (i.e. the BL does ratings of a trial right after feedback while participants do it at the beginning of a trial))
      # for the other ratings, we make a new column with otherpre-drift values for all pre-drift trials, and same for post-drift trials
      # NOTE: the first post-drift trial (phase_trial=12) still gets pre-drift value (mu & uncertainty) bc this was the prior for the post-drift learner (although it was only the mu prior, uncertainty prior was 30)
      tmp = dataBL %>% 
        group_by(ID2,sample,ID,game_number,phase) %>% mutate(across(starts_with('BL_'), ~ ifelse(phase_trial==1,NA,lag(.x)))) %>%  # shift all BL columns by one trial/row, replace first trial with NA
        ungroup() %>% mutate(
          # get mean values
          BL_ctrl_mu  = ifelse(phase_trial==1,rate_ctrl, BL_ctrl_mu), # replace mu of first phase trial by prior (first rating)
          BL_other_mu = ifelse(phase_trial==1,rate_other,ifelse(player_mu_drift==0 | phase_trial==12,BL_otherpre_mu, BL_otherpost_mu)),
          BL_self_mu  = ifelse(phase_trial==1 & phase==1,rate_self, BL_self_mu),
          # for the prior upper and lower uncertainty, I take the 90% confidence interval around the prior mean assuming the respective standard deviation (e.g. rate_self - 5 * qnorm(0.95)) => is 0.95 because we do this two-tailed for upper and lower CI
          BL_self_lo  = ifelse(phase_trial==1 & phase==1,BL_self_mu - qnorm(0.95)*sim_results$sd2fit, BL_self_lo), # sim_results$sd2fit was the prior SD i.e. initial uncertainty, qnorm gets z value
          BL_self_up  = ifelse(phase_trial==1 & phase==1,BL_self_mu + qnorm(0.95)*sim_results$sd2fit, BL_self_up),
          BL_other_lo = ifelse(phase_trial==1,BL_other_mu - qnorm(0.95)*sim_results$sd2fit, ifelse(player_mu_drift==0 | phase_trial==12,BL_otherpre_lo,BL_otherpost_lo)), 
          BL_other_up = ifelse(phase_trial==1,BL_other_mu + qnorm(0.95)*sim_results$sd2fit, ifelse(player_mu_drift==0 | phase_trial==12,BL_otherpre_up,BL_otherpost_up)),
          BL_ctrl_lo  = ifelse(phase_trial==1,BL_ctrl_mu - qnorm(0.95)*(sim_results$sd2fit/100), BL_ctrl_lo), 
          BL_ctrl_up  = ifelse(phase_trial==1,BL_ctrl_mu + qnorm(0.95)*(sim_results$sd2fit/100), BL_ctrl_up),
          # uncertainty (compute the uncertainty as the abs difference between upper and lower CI)
          BL_self_unc  = abs(BL_self_up - BL_self_lo), BL_other_unc = abs(BL_other_up - BL_other_lo), BL_ctrl_unc = abs(BL_ctrl_up - BL_ctrl_lo),
          # copy active_trial and feedback column over because I modified these for specific learners (Il and PL)
          BL_active_trial = active_trial, BL_feedback = feedback) %>%
        rowwise() %>% mutate(BL_avg_unc = mean(c(BL_self_unc,BL_other_unc,100*BL_ctrl_unc),na.rm=T),
                             BL_max_unc = ifelse(phase!=2,max(c(BL_self_unc,BL_other_unc,100*BL_ctrl_unc),na.rm=T),NA)) %>% 
        select(ID2,sample,ID,game_number,phase,phase_trial,starts_with('BL_')) %>%
        rename_with(., ~ gsub('BL_',paste0(learner2load,'_'),.x,ignore.case=F)) # now replace 'BL_' prefix in column names by the respective abbreviation of our current learner
      
      # overwrite dataP with dataP2
      dataP2 = dataP %>% left_join(tmp,by = c("game_number", "phase", "phase_trial", "ID"))
      
      # rename variables according to sample
      assign(paste0('dataB.',ifolder),dataP2)
      assign(paste0('dataBL.',ifolder),dataBL)
      
    }
  }
  
  # do some data transformations that take care of problems that are because PL doesn't have all trials
  # - take care of some duplicate columns that have odd values because the PL doesn't have all trials
  # - PL doesn't learn from active trials -> I just copy over the previous trials' PL value for those trials
  PL_vars2copy = c('PL_self_up','PL_self_lo','PL_self_mu','PL_self_unc','PL_other_up','PL_other_lo','PL_other_mu','PL_other_unc','PL_ctrl_up','PL_ctrl_lo','PL_ctrl_mu','PL_ctrl_unc','PL_avg_unc','PL_max_unc')
  dataB.mri = dataB.mri %>% select(-ID2.y,-ID2,-sample,-sample.y) %>% rename(ID2 = ID2.x, sample = sample.x) %>%
    group_by(ID,game_number,phase) %>% 
    mutate(PL_active_1after=lag(PL_active_trial), # just a helper variable for below
           across(all_of(PL_vars2copy), ~ ifelse(active_trial==1 & is.na(PL_active_trial),lag(.x),.x))) # current trial is active
  # check if we still have NA rows & need to copy PL rows over again
  stillNAs = c(dataB.mri %>% filter(phase==1) %>%pull(PL_self_up),dataB.mri %>% filter(phase==3) %>% pull(PL_ctrl_up)) %>% is.na() %>% any()
  while(stillNAs==T) {
    dataB.mri = dataB.mri %>% group_by(ID,game_number,phase) %>% 
      mutate(across(all_of(PL_vars2copy), ~ ifelse(active_trial==1 & is.na(PL_active_trial) & active_1after==1 & is.na(PL_active_1after),lag(.x),.x)))
    stillNAs = c(dataB.mri %>% filter(phase==1) %>%pull(PL_self_up),dataB.mri %>% filter(phase==3) %>%pull(PL_ctrl_up)) %>% is.na() %>% any()
  }
  
  # now do same data transformations for dataB.online
  dataB.online = dataB.online %>% select(-ID2.y,-ID2,-sample,-sample.y) %>% rename(ID2 = ID2.x, sample = sample.x) %>%
    group_by(ID,game_number,phase) %>% 
    mutate(PL_active_1after=lag(PL_active_trial), # just a helper variable for below
           across(all_of(PL_vars2copy), ~ ifelse(active_trial==1 & is.na(PL_active_trial),lag(.x),.x)))# current trial is active
  # check if we still have NA rows & need to copy PL rows over again
  stillNAs = c(dataB.online %>% filter(phase==1) %>%pull(PL_self_up),dataB.online %>% filter(phase==3) %>%pull(PL_ctrl_up)) %>% is.na() %>% any()
  while(stillNAs==T) {
    dataB.online = dataB.online %>% group_by(ID,game_number,phase) %>% 
      mutate(across(all_of(PL_vars2copy), ~ ifelse(active_trial==1 & is.na(PL_active_trial) & active_1after==1 & is.na(PL_active_1after),lag(.x),.x)))
    stillNAs = c(dataB.online %>% filter(phase==1) %>%pull(PL_self_up),dataB.online %>% filter(phase==3) %>%pull(PL_ctrl_up)) %>% is.na() %>% any()
  }

  # put together output & return the files
  out = list(dataB.mri = dataB.mri, dataBL.mri = dataBL.mri,dataB.online = dataB.online, dataBL.online = dataBL.online)
  return(out)
}