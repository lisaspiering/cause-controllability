load_data4regs <- function(path,subs2filter,reprocess=F) {
  # load data for regressions
  # path: path to the data
  # subs2filter: list of subjects to include/exclude
  # reprocess: if TRUE, reprocess the data
  
  # reload the data depending on whether we've already processed it before
  if (reprocess) {
    
    # A) load MRI behaviour
    load(file = file.path(path$data,'preprocessed_behaviour_mri.RData'))
    dataB.mri = dataB %>% filter(ID %in% subs2filter$mri.incl); dataQ.mri = dataQ %>% filter(ID %in% subs2filter$mri.incl)
    
    # B) load online data
    load(file.path(path$data,'preprocessed_behaviour_online.RData'))
    dataQ.online = dataQ.online %>% filter(!(ID %in% subs2filter$online.excl)) %>%
      # transform debrief to numerical for easier stats
      mutate(DQ.ATSelf.num  = ifelse(DQ.ATSelf  =='Very much',6, ifelse(DQ.ATSelf  =='Not at all',1,DQ.ATSelf)) %>% as.numeric(),  
             DQ.ATOther.num = ifelse(DQ.ATOther =='Very much',6, ifelse(DQ.ATOther =='Not at all',1,DQ.ATOther)) %>% as.numeric(),
             DQ.ATCtrl.num  = ifelse(DQ.ATCtrl  =='Very much',6, ifelse(DQ.ATCtrl  =='Not at all',1,DQ.ATCtrl)) %>% as.numeric())
    dataB.online = dataB.online %>% filter(!(ID %in% subs2filter$online.excl)) %>%
      mutate(
        feedback_mu        = self_mu * control_level + player_mu * (1-control_level),
        rate_upd_sumAbs    = ifelse(phase==1,abs(rate_self_updAbs + rate_other_updAbs),
                                    ifelse(phase==3, abs(rate_self_updAbs + rate_other_updAbs + 100*rate_ctrl_updAbs),NA)),
        # how much will you update on the next trial?
        rate_self_updAbsN  = ifelse(phase_trial==16,NA,lead(rate_self_updAbs)),  rate_self_updN     = ifelse(phase_trial==16,NA,lead(rate_self_upd)),
        rate_other_updAbsN = lead(rate_other_updAbs),                            rate_other_updN    = lead(rate_other_upd),
        rate_ctrl_updAbsN  = lead(rate_ctrl_updAbs),                             rate_ctrl_updN     = lead(rate_ctrl_upd),
        # errors
        rate_self_err      = abs(rate_self - self_mu),                           rate_other_err     = abs(rate_other - player_mu),
        rate_ctrl_err      = abs(rate_ctrl - control_level),                     rate_total_err     = abs(feedback_mu - rate_total))
    
    # C) load Bayesian learner results (that was fit before)
    tmp  = load_learners(dataB.mri,dataB.online,path)
    
    # D) transform data
    tmp2 = lapply(list(dataB.mri = tmp$dataB.mri,dataB.online = tmp$dataB.online), 
                  function(.x) {
                    .x %>%
                      mutate(
                        active_next                  = ifelse(phase_trial==16,NA,lead(active_trial)),
                        game_numberF                 = factor(game_number,levels = c('1','2','3','4'),ordered = T),
                        rate_sum_err                 = ifelse(phase==1,rate_self_err+rate_other_err,ifelse(phase==3,rate_self_err+rate_other_err+100*rate_ctrl_err,NA)),
                        gameS = ifelse(game_name=='slingshot', 1, 0), gameR = ifelse(game_name=='rockslide', 1, 0), gameB = ifelse(game_name=='bouncingball', 1, 0), gameL = ifelse(game_name=='lightGame', 1, 0), 
                        # compute the learner uncertainty as the abs difference between upper and lower CI
                        BL_self_unc  = abs(BL_self_up - BL_self_lo), BL_other_unc = abs(BL_other_up - BL_other_lo), BL_ctrl_unc = abs(BL_ctrl_up - BL_ctrl_lo)
                      ) %>% 
                      group_by(ID,game_number,phase) %>% mutate(
                        rate_ctrl_updAbs_prev        = lag(rate_ctrl_updAbs),
                        feedback_prev                = lag(feedback),
                        total_PE_signed_prev         = lag(total_PE_signed),
                        # rating error changes
                        rate_self_err_upd            = c(NA,diff(rate_self_err)),   rate_other_err_upd  = c(NA,diff(rate_other_err)),    rate_ctrl_err_upd  = c(NA,diff(rate_ctrl_err)),    rate_sum_err_upd= c(NA,diff(rate_sum_err)),
                        # previous ratings
                        rate_self_prev               = lag(rate_self),      rate_other_prev      = lag(rate_other),   rate_ctrl_prev  = lag(rate_ctrl),
                        # next ratings
                        rate_self_next               = lead(rate_self),     rate_other_next      = lead(rate_other),  rate_ctrl_next  = lead(rate_ctrl),
                        
                        ## LEARNER-SPECIFIC
                        BL_ctrl_unc_updAbs           = c(NA,abs(diff(BL_ctrl_unc))), IL_ctrl_unc_updAbs = c(NA,abs(diff(IL_ctrl_unc))), PL_ctrl_unc_updAbs = c(NA,abs(diff(PL_ctrl_unc))),
                        # previous uncertainties
                        BL_self_unc_prev   = lag(BL_self_unc),      IL_self_unc_prev   = lag(IL_self_unc),      PL_self_unc_prev       = lag(PL_self_unc),      
                        BL_other_unc_prev  = lag(BL_other_unc),     IL_other_unc_prev  = lag(IL_other_unc),     PL_other_unc_prev      = lag(PL_other_unc),
                        BL_ctrl_unc_prev   = lag(BL_ctrl_unc),      IL_ctrl_unc_prev   = lag(IL_ctrl_unc),      PL_ctrl_unc_prev       = lag(PL_ctrl_unc),      
                        BL_avg_unc_prev    = lag(BL_avg_unc),       IL_avg_unc_prev    = lag(IL_avg_unc),       PL_avg_unc_prev        = lag(PL_avg_unc), 
                        BL_max_unc_prev    = lag(BL_max_unc),       IL_max_unc_prev    = lag(IL_max_unc),       PL_max_unc_prev        = lag(PL_max_unc),
                        # uncertainty update change from previous to current trial; updN is from current to next trial
                        BL_self_unc_upd              = ifelse(phase_trial<2,NA,c(NA,diff(BL_self_unc))),  BL_self_unc_updN  = lead(BL_self_unc_upd),
                        BL_other_unc_upd             = ifelse(phase_trial<2,NA,c(NA,diff(BL_other_unc))), BL_other_unc_updN = lead(BL_other_unc_upd),
                        BL_ctrl_unc_upd              = ifelse(phase_trial<2,NA,c(NA,diff(BL_ctrl_unc))),  BL_ctrl_unc_updN  = lead(BL_ctrl_unc_upd),
                        BL_avg_unc_upd               = ifelse(phase_trial<2,NA,c(NA,diff(BL_avg_unc))),   BL_avg_unc_updN   = lead(BL_avg_unc_upd),
                        # BL errors
                        BL_self_mu_err               = abs(BL_self_mu - self_mu),   BL_other_mu_err              = abs(BL_other_mu - player_mu),               BL_ctrl_mu_err               = abs(BL_ctrl_mu - control_level), 
                        BL_mu_sum_err                = ifelse(phase==1,BL_self_mu_err+BL_other_mu_err,ifelse(phase==3,BL_other_mu_err+100*BL_ctrl_mu_err,NA)),
                        IL_self_mu_err               = abs(IL_self_mu - self_mu),   IL_other_mu_err              = abs(IL_other_mu - player_mu),               IL_ctrl_mu_err               = abs(IL_ctrl_mu - control_level), 
                        IL_mu_sum_err                = ifelse(phase==1,IL_self_mu_err+IL_other_mu_err,ifelse(phase==3,IL_other_mu_err+100*IL_ctrl_mu_err,NA)),
                        PL_self_mu_err               = abs(PL_self_mu - self_mu),   PL_other_mu_err              = abs(PL_other_mu - player_mu),               PL_ctrl_mu_err               = abs(PL_ctrl_mu - control_level), 
                        PL_mu_sum_err                = ifelse(phase==1,PL_self_mu_err+PL_other_mu_err,ifelse(phase==3,PL_other_mu_err+100*PL_ctrl_mu_err,NA)),
                        # BL rating updates
                        BL_self_mu_upd               = c(NA,diff(BL_self_mu)),
                        BL_other_mu_upd              = c(NA,diff(BL_other_mu)),
                        BL_ctrl_mu_upd               = c(NA,diff(BL_ctrl_mu)),
                        # BL total PEs
                        BL_rate_total           = ifelse(phase==1, BL_self_mu * control_level + BL_other_mu * (1-control_level), 
                                                         ifelse(phase==3,self_mu * BL_ctrl_mu + BL_other_mu * (1-BL_ctrl_mu), NA)),
                        BL_total_PE_signed      = feedback - BL_rate_total,
                        BL_total_PE_signed_prev = ifelse(phase!=2,lag(BL_total_PE_signed),NA)
                      ) %>% ungroup() %>%
                      mutate(across(c(active_trial,active_1after,active_2after,active_3after,active_next,player_mu_drift),as.factor))
                  })
    
    dataB.mri     = tmp2$dataB.mri;   
    dataB.online  = tmp2$dataB.online;
    
    
    # put everything together and return
    out = list(dataB.mri = dataB.mri,dataB.online = dataB.online, dataQ.mri = dataQ.mri,dataQ.online = dataQ.online)
    save(out,file = file.path(path$gitlab,'behavioural-data','preprocessed_data_4regs.RData'))
  } else {
    load(file = file.path(path$gitlab,'behavioural-data','preprocessed_data_4regs.RData'))
  }
  return(out)
}