fit_sub_behaviour = function(fit.data=fit.data,fit.vars=fit.vars,fit.model=fit.model,fit.name=fit.name,path=path,my.priors=my.priors,max.tries=5){
  # PART 1
  # inputs: fit.data - data to be fitted (contains all subjects)
  #         fit.model - stan model, already compiled
  #         path - path where the subject folders are saved (i.e. it wants path$save variable)
  #         fit.name - file name name under which each subject's regressor estimates should be saved
  # outputs: pars - parameter estimate for the person that was fitted => this is saved in a file (i.e. not an output variable of calling the function directly)
  # Lisa Spiering, 2025
  
  data.sample = ifelse(grepl('.online',fit.name), 'online', 'mri') # find out which data sample we're analysing so that we can save the brms output in the right subfolder 
  subfolder   = paste0('indivs_',data.sample) 
  
  suppressWarnings(rm(dfT))
  foreach(is=fit.data$ID %>% unique() %>% sort()) %dopar% { # loop over participants
    dfT = fit.data %>% dplyr::filter(ID==is) %>% mutate_at(fit.vars$toscale,scale)
    
    # make a new subject folder if not yes exists
    ipath = file.path(path$save,subfolder,paste0('s',str_pad(is, 2, pad = "0")))
    if (!file.exists(ipath)) {dir.create(ipath)}
    
    if (nrow(dfT)>0){
      count_divergent = 1 #initialize parameters to keep track of fitting to criterion
      my.rhat         = 5
      my.adapt_delta  = 0.8
      my.iter         = 8000
      my.cores        = 1 # quicker for lowest number of iterations & adapt_delta to run it on 1 core only, we parallelise as soon as we need to rerun this with more iterations
      try.count       =  1
      suppressWarnings(rm(out)) # rm output if it exists, and suppress warning if it doesn't exist yet
      while (count_divergent >0 || my.rhat>=1.1){
        out = update(fit.model,newdata=dfT,iter=my.iter,control=list(adapt_delta=my.adapt_delta),cores=my.cores, chains=4,refresh=0, prior = my.priors) # 
        # check divergences, rhat and rerun if necessary
        np              = nuts_params(out)
        count_divergent = sum(subset(np,Parameter=="divergent__")$Value)
        my.rhat         = rhat(out)[names(rhat(out))!=c('lp__')] %>% max() # don't check rhat for lp__
        my.iter         = my.iter*1.5
        my.adapt_delta  = my.adapt_delta + (1-my.adapt_delta)*0.5
        my.cores        = 1 # reset 4 to 1 bc I do parallel computing above
        try.count       = try.count+1
        if (try.count>max.tries){
          break
        }
      }
      # Note down whether fitting was successful and extract parameters if successful
      par.names = rownames(fixef(out)) # extract parameter names so that the code can be used for each model
      
      if (try.count<=max.tries){
        pars = data.frame(Estimate=fixef(out)[,1]) %>% t() %>% as.data.frame() %>% # extract fixed effects / 1 is the Estimate (other columns are the est. error and 2.5% and 97.5% interval)
          mutate(ID = unique(dfT$ID), fitted.ok.flag = 1)
      } else {
        pars           = data.frame(Estimate=fixef(out)[,1]) %>% t() %>% as.data.frame() # first extract fixed effects as well but just to get the column names, then replace the ill-fitted parameters with NA
        pars[1,]       = NA # then replace the estimated beta values by NA because the model fit didn't work for this subject
        pars           = pars %>% mutate(ID = unique(dfT$ID), fitted.ok.flag = 0)
      }
      
      save(pars,file=file.path(ipath, paste0(fit.name,'.RData'))) # str_pad creates a string such as s0001
    }
  }
  
  # PART 2
  # this part loads all the subject regression fits, and then puts them together into one dataframe with one row per subject, and columns for each regressor. The dataframe is then saved in the all_subjects folder in path$save
  # inputs: fit.data = data frame with all the data (e.g. dataB), fit.name = filenames of the subject files which we load, path = path variable where the subject files are stored
  # output: df.all.subjects = dataframe with columns for all regressors estimates, and rows for each subject
  
  df.all.subjects = data.frame()
  for (is in unique(fit.data$ID)){
    thisFile = file.path(path$save,subfolder,paste0('s',str_pad(is, 2, pad = "0")),paste0(fit.name,'.RData'))
    if (file.exists(thisFile)){
      load(thisFile)
      df.all.subjects = rbind(df.all.subjects,as.data.frame(pars))
    }
  }
  
  # PART 3
  # do t.test of individual regressors against 0
  subs.notfittedok = df.all.subjects %>% filter(fitted.ok.flag==0) %>% pull(ID)
  ttests = lapply(df.all.subjects %>% filter(!(ID%in%subs.notfittedok)) %>% select(-ID,-fitted.ok.flag),t.test,rep(0,nrow(df.all.subjects%>%filter(!(ID%in%subs.notfittedok)))))
  ttests$exclIDnotfittedok = subs.notfittedok
  
  # put all the results together
  results = list(df = df.all.subjects, ttests = ttests, model = fit.model) # put results in a list
  assign(fit.name,results) # rename results variable
  save(list = fit.name,file=file.path(path$save,subfolder,'all_subs',paste0(fit.name,'.RData'))) # list = fit.name means we save the variable that has the name of the string saved under fit.name (i.e. we save e.g. fit_self_upAbs... rather than the variable fit.name)
  
  
  return(results)
}