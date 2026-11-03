data{
  int<lower=1> N;
  int<lower=1> N_pred;
  vector[N] swb;
  // vector[N] temp;
  // vector[N] avgba;
  int abun[N];
  
  real min_swb;
  real max_swb;
  real pres_swb_dat;
  
  
  vector[N_pred] swb_pred;

}

parameters{
  real a;
  real b_swb;
  real b2_swb;
  // real b_temp;
  // real b2_temp;
  // real b_ba;
  
  real<lower=0> phi;
  
}

transformed parameters{
  vector[N] eta;
    
  for(i in 1:N){
   eta[i] = a + b_swb*swb[i] + b2_swb*square(swb[i]);
      // + b_temp*temp[i] + b2_temp*square(temp[i]) + b_ba*avgba[i];
  }
}

model{
  a ~ normal(0, 10);
  
  b_swb ~ normal(0, 10);
  b2_swb ~ normal(0, 10);
  // b_temp ~ normal(0, 10);
  // b2_temp ~ normal(0, 10);
  // b_ba ~ normal(0, 10);
  
  phi ~ gamma(0.01, 0.01);
  
  abun ~ neg_binomial_2_log(eta, phi);
  
}

generated quantities{
  // int y_rep[N];
  real eta_pred[N_pred];
  
  real wet_swb;
  real dry_swb;
  real pres_swb;
  real par_swb;
  real opt_swb;
  real n_hat_wet;
  real n_hat_dry;
  real n_hat_pres;
  real n_hat_opt;
  
  // real log_lik[N];
  // for(i in 1:N) log_lik[i] = neg_binomial_2_log_lpmf(abun[i] | eta[i], phi);


  for(i in 1:N_pred){
   eta_pred[i] = exp(a + b_swb*swb_pred[i] + b2_swb*square(swb_pred[i]));
  }

  dry_swb = min_swb;
  wet_swb = max_swb;
  pres_swb = pres_swb_dat;
  par_swb = -b_swb/(2*b2_swb);
    if ((par_swb <= wet_swb) && (par_swb >= dry_swb)){
      opt_swb = par_swb;
    }  else if (par_swb >= wet_swb){
        opt_swb = wet_swb;
      }  else {
        opt_swb = dry_swb;
      }

  n_hat_wet = exp(a + b_swb*wet_swb + b2_swb*square(wet_swb));
  n_hat_dry = exp(a + b_swb*dry_swb + b2_swb*square(dry_swb));
  n_hat_pres = exp(a + b_swb*pres_swb + b2_swb*square(pres_swb));
  n_hat_opt = exp(a + b_swb*opt_swb + b2_swb*square(opt_swb));

}

