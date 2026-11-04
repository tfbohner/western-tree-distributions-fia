data{
  int<lower=1> N;
  int<lower=1> N_pred;
  vector[N] swb;
  int pres[N];
  
  real min_swb;
  real max_swb;
  
  vector[N_pred] swb_pred;
}

parameters{
  real a;
  real b_swb;
  real b2_swb;
  
}

transformed parameters{
  vector[N] theta; 
  for(i in 1:N){
   theta[i] = a + b_swb*swb[i] + b2_swb*square(swb[i]);
  }
}

model{
  a ~ normal(0, 10);
  b_swb ~ normal(0, 10);
  b2_swb ~ normal(0, 10);

  pres ~ bernoulli_logit(theta);
}

generated quantities{
  real theta_pred[N_pred];
  
  real par_swb;
  real opt_swb;
  real thetamin_swb;
  real thetamax_swb;
  real thetaopt_swb;
  
  

  
  
    for(i in 1:N_pred){
   theta_pred[i] = a + b_swb*swb_pred[i] + b2_swb*square(swb_pred[i]);
  }

par_swb = -b_swb/(2*b2_swb);
    if ((min_swb <= par_swb) && (par_swb <= max_swb)){
      opt_swb = par_swb;
    }  else if (min_swb >= par_swb){
        opt_swb = min_swb;
      }  else {
        opt_swb = max_swb;
      }
    thetamin_swb = a + b_swb*min_swb + b2_swb*square(min_swb);
    thetamax_swb = a + b_swb*max_swb + b2_swb*square(max_swb);
    thetaopt_swb = a + b_swb*opt_swb + b2_swb*square(opt_swb);
 
}

