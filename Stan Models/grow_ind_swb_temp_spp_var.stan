data {
  int<lower=1> N; 
  int<lower=1> N_pred;
  int<lower=1> N_plts;
  int<lower=1> plot[N];
  vector[N] dia;
  vector[N] swb;
  vector[N] avgba;
  vector[N] adcount;
  vector[N] temp;
  vector[N] y; 

  vector[N_pred] swb_pred;
  real min_swb;
  real max_swb;
  real pres_swb_dat;
  // real abun_swb_dat;
}

parameters {
  real a[N_plts];  
  real a_sp;  
  real b_dia;
  real b_swb;
  real b_ba;
  real b_adcount;
  real b_temp;
  real b2_dia;
  real b2_swb;
  real b2_temp;
  
  real sigma_dia; 
  real<lower=0.00000001> sigma_c; 
}

transformed parameters {
  vector[N] mu;
  vector[N] sigma; 

  for(i in 1:N){
    mu[i] = a[plot[i]] + b_dia*dia[i] + b_swb*swb[i] + b_ba*avgba[i] + b_adcount*adcount[i] + b_temp*temp[i] + b2_dia*square(dia[i]) + b2_swb*square(swb[i]) + b2_temp*square(temp[i]); 
  }
  
  for(i in 1:N){
    sigma[i] = exp(sigma_dia*dia[i]); 
  }
}

model {
  a ~ normal(a_sp, sigma_c);
  a_sp ~ normal(0, 10);

  b_dia ~ normal(0, 10); 
  b_swb ~ normal(0, 10); 
  b_temp ~ normal(0, 10); 
  b_ba ~ normal(0, 10); 
  b_adcount ~ normal(0, 10); 
  b2_dia ~ normal(0, 10);
  b2_swb ~ normal(0, 10);
  b2_temp ~ normal(0, 10);
  
  sigma_c ~ normal(0, 1);
  sigma_dia ~ normal(0, 10);

  //Likelihood
  y ~ normal(mu, sigma); 
}

generated quantities{
  real wet_swb;
  real dry_swb;
  real pres_swb;
  // real abun_swb;
  real par_swb;
  real opt_swb;
  real g_hat_wet;
  real g_hat_dry;
  real g_hat_pres;
  real g_hat_abun;
  real g_hat_opt;

  real mu_pred[N_pred];
  real log_lik[N];
  
  for(i in 1:N) log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);

  for(i in 1:N_pred){
       mu_pred[i] = a_sp + b_swb*swb_pred[i] + b2_swb*square(swb_pred[i]);
       }

    dry_swb = min_swb;
    wet_swb = max_swb;
    pres_swb = pres_swb_dat;
    // abun_swb = abun_swb_dat;
    par_swb = -b_swb/(2*b2_swb);
    if ((par_swb <= wet_swb) && (par_swb >= dry_swb)){
      opt_swb = par_swb;
    }  else if (par_swb >= wet_swb){
        opt_swb = wet_swb;
      }  else {
        opt_swb = dry_swb;
      }
    g_hat_wet = a_sp + b_swb*wet_swb + b2_swb*square(wet_swb);
    g_hat_dry = a_sp + b_swb*dry_swb + b2_swb*square(dry_swb);
    g_hat_pres = a_sp + b_swb*pres_swb + b2_swb*square(pres_swb);
    // g_hat_abun = a_sp + b_swb*abun_swb + b2_swb*square(abun_swb);
    g_hat_opt = a_sp + b_swb*opt_swb + b2_swb*square(opt_swb);

}
