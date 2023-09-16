
data {

	int K; /* Number of studies */
	int J; /* sensitivity and specificity = 2*/
	int TP[K];
	int FP[K];
	int FN[K];
	int TN[K];
	int prior_type;

}

transformed data {
	int Y_pos[K];
	int Y_neg[K];


	for(i in 1:K){
		Y_pos[i] = TP[i] + FN[i]; // num positive
		Y_neg[i] = TN[i] + FP[i]; // num negative
		}
}

parameters{
	vector[J] mu;
	vector[J] theta[K];
	vector<lower=0>[J] sigma; // between-study sd
	cholesky_factor_corr[J] Lcorr;


}

transformed parameters{
	vector[J] P[K] = inv_logit(theta);
	
}

model {
	
	for (i in 1:K){
		TP[i] ~ binomial(Y_pos[i], P[i,1]);
		TN[i] ~ binomial(Y_neg[i], P[i,2]);
	}


	theta ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, Lcorr));
	

	// priors
    if (prior_type == 1) {
     	mu ~ normal(0,2);
	    sigma ~ uniform(0, 5);
	    Lcorr ~ lkj_corr_cholesky(4);
    } else if (prior_type == 2) {
     	mu ~ normal(0,5);
	    sigma ~ uniform(0, 10);
	    Lcorr ~ lkj_corr_cholesky(6);
    } 
	
}

generated quantities {
	vector<lower=0, upper=1>[J] pool_P;
	matrix[J,J] Omega_b;
	vector[2] theta_pred;
	vector[2] sesp_pred;
  vector[K] log_lik;

  Omega_b  = multiply_lower_tri_self_transpose(Lcorr); // Correlation matrix

	pool_P = inv_logit(mu); //posterior pooled sensitivity and specificity

	// predictive distributions (new)
	theta_pred = multi_normal_cholesky_rng(mu, diag_pre_multiply(sigma, Lcorr));
	sesp_pred = inv_logit(theta_pred);
	
	// For DIC
	for (i in 1:K) {
    	log_lik[i] = binomial_lpmf(TP[i] | Y_pos[i], P[i,1]) +
                 	 binomial_lpmf(TN[i] | Y_neg[i], P[i,2]);
  }
}
