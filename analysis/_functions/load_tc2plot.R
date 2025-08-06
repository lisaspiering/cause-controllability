load_tc2plot <- function(regnames,roi2plot,allBeta) {
  # convenience plotting variables
  # Lisa Spiering, 2025
  time = seq(from=-2,to=15, length.out=250) # 250 samples, -2 to 15 is our plotted time epoch
  
  betas2plot = list()
  for (eoi2plot in names(regnames)) {
    # name regressors in allBeta
    tc2plot = allBeta[[roi2plot]][1:length(regnames)][[eoi2plot]][[1]]
    
    # put participants from being in separate slices in 3rd dimension into one long 2D data frame
    tc_all_subs = data.frame()
    for (s in c(1:dim(tc2plot)[3])) {
      sbeta = tc2plot[,,s] %>% t() %>% as.data.frame() %>% set_colnames(regnames[[eoi2plot]]) %>% mutate(ID = subs2filter$mri.incl[s], time = time)
      tc_all_subs = rbind(tc_all_subs,sbeta)
    }
    tc_all_subs = tc_all_subs %>% pivot_longer(cols=-c(ID,time),names_to='var',values_to='beta')
    
    # get group stats
    tc_mu_se = tc_all_subs %>% group_by(time,var) %>% dplyr::summarise(beta_mu = mean(beta), beta_se = sd(beta) / sqrt(dim(tc2plot)[3]),
                                                                       beta_seL = beta_mu - beta_se, beta_seH = beta_mu + beta_se) %>% ungroup()
    
    # store in list
    betas2plot[[eoi2plot]] = tc_mu_se
  }
  return(betas2plot)
}