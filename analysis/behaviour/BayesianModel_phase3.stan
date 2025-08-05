data {
  int<lower=0> nrow;
  real feedback[nrow];
  int active_trial[nrow];
  int player_mu_drift[nrow];
  real self_mu[nrow]; // here we tell BL the true self (phase 3)
  real self_sd[nrow];
  real PriorOtherMean; 
  real PriorOtherSD; 
  real PriorControlMean;
  real PriorControlSD;
  real PriorOtherMeanDrift; 
}

parameters {
  real<lower=0,upper=100> RateOtherMean;
  real<lower=0> RateOtherSD;
  real<lower=0,upper=1> RateControlMean;
  real<lower=0,upper=1> RateControlSD;
  real<lower=0,upper=100> RateOtherMeanDrift;
  real<lower=0> RateOtherSDDrift;
}

model {
  // specify priors
  RateOtherMean  ~ normal(PriorOtherMean, PriorOtherSD);
  RateOtherSD    ~ normal(0,5);
  RateControlMean ~ normal(PriorControlMean, PriorControlSD);
  RateOtherMeanDrift ~ normal(PriorOtherMeanDrift, PriorOtherSD); //prior for other player post-drift
  RateOtherSDDrift   ~ normal(0,5);
  
  // specify model
  for (irow in 1:nrow) { // go through all the past trials
  // phase 3: in this model, because we don't run the previous two phases, we just assume that the self level in phase 3 is 'known'
  
  if (player_mu_drift[irow]==0) {
    if (active_trial[irow]==1) {
      feedback[irow] ~ normal((1-RateControlMean) * RateOtherMean, sqrt( ((1-RateControlMean)^2) * (RateOtherSD^2)));
    } else {
      feedback[irow] ~ normal(RateControlMean * self_mu[irow] + (1-RateControlMean) * RateOtherMean,sqrt( (RateControlMean^2) * (self_sd[irow]^2) + ((1-RateControlMean)^2) * (RateOtherSD^2)));
    }
  } else {
    if (active_trial[irow]==1) {
      feedback[irow] ~ normal((1-RateControlMean) * RateOtherMeanDrift, sqrt( ((1-RateControlMean)^2) * (RateOtherSDDrift^2)));
    } else {
      feedback[irow] ~ normal(RateControlMean * self_mu[irow] + (1-RateControlMean) * RateOtherMeanDrift,sqrt( (RateControlMean^2) * (self_sd[irow]^2) + ((1-RateControlMean)^2) * (RateOtherSDDrift^2)));
    }
  }
  }
}
