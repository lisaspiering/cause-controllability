transform_data_l1 <- function(dataB.mri,path,subs2filter) {
  # this function does some data transformations specific to the fmri analysis
  # we mainly select the columns of dataB.mri that we're interested in, filter out Self-Only trials (phase 2), and then recompute some onsets, durations, and RTs
  # input: dataB.mri
  # output: dataB.mri
  # Lisa Spiering, 2025
  dataB.mri = out$dataB.mri %>% filter(phase!=2) %>%
    select(ID,game_number,phase,game_trial,phase_trial,gameS,gameB,gameR,gameL,starts_with('rate_'),starts_with('active_'),self_feedback,feedback,feedback_prev,starts_with('total_PE'),player_mu,self_mu,control_level,player_mu_drift,starts_with('timing_'),starts_with('onset'),starts_with('offset'),
           contains('BL_avg_unc'),contains('BL_max_unc'), contains('BL_self_unc'),contains('BL_other_unc'),contains('BL_ctrl_unc'),-contains('curprev'),-contains('AbsN2'),perferror_abs) %>%
    mutate(constant = 1, # constant column for all roi designs
           across(starts_with('active_'), ~ as.integer(.x)-1), # otherwise entries are saved as 1 or 2 for matlab for some reason
           # onsets
           onset_Dc1  = timing_getRating,       offset_Dc1  = timing_ratedself-1, # for offset of control, take 1 sec after looking at ctrl
           onset_D13  = timing_getRating,       offset_D13  = timing_confirmRating,
           onset_Ds1  = timing_ratedself-1   ,  offset_Ds1  = timing_ratedself,        # we substract one second (1000ms) from when they finished rating self for the onset
           onset_Do13 = timing_ratedself,       offset_Do13 = timing_ratedother, 
           onset_Dc3  = timing_getRating,       offset_Dc3  = timing_ratedcontrol,
           onset_Ds3  = timing_ratedcontrol,    offset_Ds3  = timing_ratedself,
           onset_A13  = timing_playMinigame,    offset_A13  = timing_finishedMinigame,
           onset_O13  = timing_giveFeedback,    offset_O13  = timing_confirmFeedback,
           # durations
           dur_D13    = abs( offset_D13  - onset_D13),
           dur_Dc1    = abs( offset_Dc1  - onset_Dc1), 
           dur_Ds1    = abs( offset_Ds1  - onset_Ds1), # is always 1s at the moment bc we don't have exact onset of self rating
           dur_Dcs1   = abs( offset_Ds1  - onset_Dc1),
           dur_Do13   = abs( offset_Do13 - onset_Do13), 
           dur_Dc3    = abs( offset_Dc3  - onset_Dc3), 
           dur_Ds3    = abs( offset_Ds3  - onset_Ds3),
           dur_A13    = abs( offset_A13  - onset_A13),
           dur_O13    = abs( offset_O13  - onset_O13), 
           # log reaction times
           RT_Do      = log(dur_Do13),  
           RT_Dcs1    = log(dur_Dcs1), 
           RT_Ds3     = log(dur_Ds3),  
           RT_Dc3     = log(dur_Dc3),
           RT_A13     = log(dur_A13),   
           RT_O13     = log(dur_O13)
    )
  
  return(dataB.mri)
}