data {
  int<lower=1> N; 
  int<lower=1> N_spp; 
  int<lower=1> spp[N];
  int sdlgs[N]; 
  int<lower=1> N_pred;
  
  real swb[N];
  real avgba[N];
  real adcount[N];
  real temp[N];

  real min_swb[N_spp];
  real max_swb[N_spp];
  real pres_swb_dat[N_spp];
  vector[N_pred] swb_pred;
}

parameters {
  real a_pres[N_spp];
  real a_abund[N_spp];
  
  real b_ba[N_spp];  
  real b_adcount[N_spp];  
  real b_temp[N_spp];
  real b2_temp[N_spp];
  real b_swb[N_spp];  
  real b2_swb[N_spp];  
  
  real<lower=0> phi[N_spp];

}

transformed parameters {
  vector[N] mu;
  vector[N] theta;
  
  for(i in 1:N){
  
  theta[i] = a_pres[spp[i]] ;
 
    mu[i] = a_abund[spp[i]] + b_adcount[spp[i]]*adcount[i] + b_ba[spp[i]]*avgba[i]+ b_swb[spp[i]]*swb[i] + b2_swb[spp[i]]*(swb[i])^2  + b_temp[spp[i]]*temp[i] + b2_temp[spp[i]]*(temp[i]^2); 

  }
  }

model {
 
  a_abund ~ normal(0,10);
  a_pres ~ logistic(0,1);
  
  b_ba ~ normal(0, 10);   
  b_adcount ~ normal(0, 10);  
  b_swb ~ normal(0, 10);  
  b2_swb ~ normal(0, 10); 
  b_temp ~ normal(0, 10);  
  b2_temp ~ normal(0, 10); 
  phi ~ gamma(0.01, 0.01);
  
      
            for (i in 1:N) {
          if (sdlgs[i] == 0) {
            target += log_sum_exp(bernoulli_logit_lpmf(1 | theta[i]), 
                               bernoulli_logit_lpmf(0 | theta[i]) + 
                               neg_binomial_2_log_lpmf(sdlgs[i] | mu[i], phi[spp[i]]));
          } else {
            target += bernoulli_logit_lpmf(0 | theta[i]) +
                               neg_binomial_2_log_lpmf(sdlgs[i] | mu[i], phi[spp[i]]);
          }
      }
  
   }

generated quantities{

  real wet_swb[N_spp];
  real dry_swb[N_spp];
  real pres_swb[N_spp];
  real par_swb[N_spp];
  real opt_swb[N_spp];
  real r_hat_wet[N_spp];
  real r_hat_dry[N_spp];
  real r_hat_pres[N_spp];
  real r_hat_opt[N_spp];

  real mu_pred[N_spp, N_pred];

  for(s in 1:N_spp){
    for(i in 1:N_pred){
      mu_pred[s,i] = a_abund[s] + b_swb[s]*swb_pred[i] + b2_swb[s]*square(swb_pred[i]);
       }}

  for(s in 1:N_spp){
    wet_swb[s] = max_swb[s];
    dry_swb[s] = min_swb[s];
    pres_swb[s] = pres_swb_dat[s];
    par_swb[s] = -b_swb[s]/(2*b2_swb[s]);
    if ((par_swb[s] <= wet_swb[s]) && (par_swb[s] >= dry_swb[s])){
      opt_swb[s] = par_swb[s];
    }  else if (par_swb[s] >= wet_swb[s]){
        opt_swb[s] = wet_swb[s];
      }  else {
        opt_swb[s] = dry_swb[s];
      }
    r_hat_wet[s] = a_abund[s] + b_swb[s]*wet_swb[s] + b2_swb[s]*square(wet_swb[s]);
    r_hat_dry[s] = a_abund[s] + b_swb[s]*dry_swb[s] + b2_swb[s]*square(dry_swb[s]);
    r_hat_pres[s] = a_abund[s] + b_swb[s]*pres_swb[s] + b2_swb[s]*square(pres_swb[s]);
    r_hat_opt[s] = a_abund[s] + b_swb[s]*opt_swb[s] + b2_swb[s]*square(opt_swb[s]);
  }

}
