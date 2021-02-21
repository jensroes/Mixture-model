/* 
Finite mixture model
This model is adapted from Vasishth et al. 2007 
Random intercepts for subj and items
Mixture on simple and complex condition
*/

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1, upper=2> condition[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

//	int<lower=1> I;                  //number of items
//	int<lower=1, upper=I> items[N];   //items id
}


parameters {
	real<lower=0> beta;			// distributions
	real<lower=0> delta;			// distribution + extra component
	real<lower=0> sigma;		// residual sd

	real<lower=0> sigma_diff;
	vector[2] theta;

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

//	vector[I] w; //items intercepts
//	real<lower=0> sigma_w;//items sd
}


transformed parameters{
  real<lower=0> mu[N];
	real<lower=0> sigmap_e;
	real<lower=0> sigma_e;
	vector[2] prob = inv_logit(theta);
	matrix[2,2] log_theta;
	matrix[N,2] log_theta_condition;

  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;

	log_theta[1, 1] = log_inv_logit(theta[1]);
  log_theta[1, 2] = log1m_inv_logit(theta[1]);
	log_theta[2, 1] = log_inv_logit(theta[2]);
  log_theta[2, 2] = log1m_inv_logit(theta[2]);


  for(n in 1:N){
    mu[n] = beta + u[subj[n]];// + w[items[n]];
    for(i in 1:2){
      log_theta_condition[n,i] = log_theta[condition[n],i]; 
    }
  }
}

model {
  vector[2] lp_parts[N];

  // Priors
	beta ~ cauchy(0, 2.5);
  delta ~ normal(0, 100);

  sigma_diff ~ normal(0,1);
  sigma ~ cauchy(0, 2.5);
  
  theta ~ normal(0, 1);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

//  sigma_w ~ normal(0,2.5);
//  w ~ normal(0, sigma_w); //items random effects


  // Likelihood	
	for(n in 1:N){
    lp_parts[n, 1] = log_theta_condition[n,1] + lognormal_lpdf(y[n] | mu[n] + delta, sigmap_e); 
    lp_parts[n, 2] = log_theta_condition[n,2] + lognormal_lpdf(y[n] | mu[n], sigma_e);
    target += log_sum_exp(lp_parts[n]); 
	}
}


generated quantities{
  vector[N] log_lik;
	vector[N] y_tilde;
  real<lower=0,upper=1> theta_long; 
  real beta2 = beta + delta;
  
  // likelihood: 
  for(n in 1:N){
    	log_lik[n] = log_sum_exp(
                  log_theta_condition[n,1] + lognormal_lpdf(y[n] | mu[n] + delta, sigmap_e), 
                  log_theta_condition[n,2] + lognormal_lpdf(y[n] | mu[n], sigma_e)
          );
    	    theta_long = bernoulli_rng(prob[condition[n]]); 
          if(theta_long) { 
              y_tilde[n] = lognormal_rng(mu[n] + delta, sigmap_e);
          }
          else{
              y_tilde[n] = lognormal_rng(mu[n], sigma_e);
          }
        }
}

