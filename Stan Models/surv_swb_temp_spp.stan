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
  int surv[N]; // 1 = survived to next census

  vector[N_pred] swb_pred;
  real min_swb;
  real max_swb;
  real pres_swb_dat;
  
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

  real<lower=0.00000001> sigma_c; // small lower bound keeps the sampler off zero

}

transformed parameters {
  vector[N] pr;

  for(i in 1:N){
    pr[i] = a[plot[i]] + b_dia*dia[i] + b_swb*swb[i] + b_ba*avgba[i] + b_adcount*adcount[i] + b_temp*temp[i] + b2_dia*square(dia[i])     + b2_swb*square(swb[i]) + b2_temp*square(temp[i]); 
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
  
  sigma_c ~ student_t(3, 0, 10);

  //Likelihood
  surv ~ bernoulli_logit(pr); 
}

generated quantities{
  real wet_swb;
  real dry_swb;
  real pres_swb;
  real par_swb;
  real opt_swb;
  real s_hat_wet;
  real s_hat_dry;
  real s_hat_pres;
  real s_hat_abun;
  real s_hat_opt;
  real pr_pred[N_pred];

  for(i in 1:N_pred){
    pr_pred[i] = a_sp + b_swb*swb_pred[i] + b2_swb*square(swb_pred[i]);
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

    s_hat_wet = a_sp + b_swb*wet_swb + b2_swb*square(wet_swb);
    s_hat_dry = a_sp + b_swb*dry_swb + b2_swb*square(dry_swb);
    s_hat_pres = a_sp + b_swb*pres_swb + b2_swb*square(pres_swb);
    s_hat_opt = a_sp + b_swb*opt_swb + b2_swb*square(opt_swb);

}
