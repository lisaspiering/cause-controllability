data {
  int<lower=0> nrow;
  real feedback[nrow];
  int active_trial[nrow];
  int player_mu_drift[nrow];
  real control_level[nrow];
  real PriorSelfMean;
  real PriorSelfSD;
  real PriorOtherMean; 
  real PriorOtherSD;
  real PriorOtherMeanDrift; 
}

parameters {
  real<lower=0,upper=100> RateSelfMean;
  real<lower=0> RateSelfSD;
  real<lower=0,upper=100> RateOtherMean;
  real<lower=0> RateOtherSD;
  real<lower=0,upper=100> RateOtherMeanDrift;
  real<lower=0> RateOtherSDDrift;
}


model {
  // specify priors
  RateSelfMean       ~ normal(PriorSelfMean, PriorSelfSD); // PriorSelfSD is initial belief uncertainty
  RateSelfSD         ~ normal(0,5); // SD of your performance
  RateOtherMean      ~ normal(PriorOtherMean, PriorOtherSD); // prior for other pre-drift
  RateOtherSD        ~ normal(0,5);
  RateOtherMeanDrift ~ normal(PriorOtherMeanDrift, PriorOtherSD); //prior for other player post-drift
  RateOtherSDDrift   ~ normal(0,5);
  
  // specify model
  for (irow in 1:nrow) { // go through all the past trials
  
  if (player_mu_drift[irow]==0) {
    if (active_trial[irow]==1) {
      feedback[irow] ~ normal((1-control_level[irow]) * RateOtherMean, sqrt( ((1-control_level[irow])^2) * (RateOtherSD^2)));
    } else {
      feedback[irow] ~ normal(control_level[irow] * RateSelfMean + (1-control_level[irow]) * RateOtherMean, sqrt( (control_level[irow]^2) * (RateSelfSD^2) + ((1-control_level[irow])^2) * (RateOtherSD^2)));
    }
  } else {
    if (active_trial[irow]==1) {
      feedback[irow] ~ normal((1-control_level[irow]) * RateOtherMeanDrift, sqrt( ((1-control_level[irow])^2) * (RateOtherSDDrift^2)));
    } else {
      feedback[irow] ~ normal(control_level[irow] * RateSelfMean + (1-control_level[irow]) * RateOtherMeanDrift, sqrt( (control_level[irow]^2) * (RateSelfSD^2) + ((1-control_level[irow])^2) * (RateOtherSDDrift^2)));
    }
  }
  
  }
}

