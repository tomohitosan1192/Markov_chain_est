#################################################################################
# STAN code to estimate Markov chain dynamics with RCS data
#
# Author: Tomohito Okabe, Hitotsubashi University, email: t-okabe@ier.hit-u.ac.jp
# For details, see Okabe, T. and Nogiwa D.,
# ``Estimation of Unobserved Dynamics of Individual Partisanship: 
# A Bayesian Approach", XXXXX, XX (20XX): - .
#################################################################################

functions{ 
  matrix matrix_pow(matrix a, int n);
  matrix matrix_pow(matrix a, int n) {
        if (n == 0)
          return diag_matrix(rep_vector(1, rows(a)));
        else
          return a *  matrix_pow(a, n - 1);
      }


}

data {
	int<lower=0> N;                     
  // Total # Observations
  int<lower=2> K;                     // Total # Responses K=3
  int<lower=1> D;                     // # of individual characteristics
  int<lower=1, upper=3> y[N];         // Response (Discrete Number/ 1,2,3)
  int<lower=0> period[N];             // Time Period Index  
  vector[D] x[N];                     // Explanatory Variables
  simplex[3] Ps_init; // Initial Distribution
}

transformed data{
  row_vector[D] zeros_beta;
  matrix<lower=0, upper=1> [D,D] A;
 
  zeros_beta = rep_row_vector(0, D);

  A = diag_matrix(rep_vector(1,D));
}

parameters {
 row_vector[D] beta_1;  
 vector[((K-1)*D)] beta_2_raw;
 row_vector[D] beta_3;  
 cov_matrix[(K-1)] Sigma_beta;

}

transformed parameters {
  matrix[K, D] beta_2; 
  cov_matrix[((K-1)*D)] var_beta_raw;

  // beta K-1 handling
  beta_2[2, :] = rep_row_vector(0, D);

  for (i in 1:(K-1)){
      for (j in 1:D){
        beta_2[(K-1)*(i-1)+1, j] = beta_2_raw[D*(i-1)+j];
      }
  }
    
  // Kronecker product
    for (i in 1:(K-1))
      for (j in 1:(K-1))
        for (k in 1:D)
          for (l in 1:D)
            var_beta_raw[D*(i-1)+k,D*(j-1)+l] = Sigma_beta[i,j]*A[k,l];
  
  }
      



model {
  vector[K] Ps[N]; 
  vector[K] phi[N,K];  // Matrix Elements     
  matrix[K,K] Mmatrix[N]; 
  int nu;
  matrix[(K-1), (K-1)] V;

  // Priors 
  nu = 2+3;
  V = nu * diag_matrix(rep_vector(1,(K-1)));

  for (d in 1:D){
    beta_1[d] ~ normal(0, 5);
    beta_3[d] ~ normal(0, 5);
  }

  Sigma_beta ~ inv_wishart(nu, V);

  beta_2_raw ~ multi_normal(rep_vector(0, ((K-1)*D) ), var_beta_raw);   

  // Likelihood
	for(n in 1:N){
         phi[n,1,3] = 0;
         phi[n,1,2] = inv_logit(beta_1 * x[n]) ;
         phi[n,1,1] = 1 -  phi[n,1,2]; 
 
         phi[n,3,1] = 0;
         phi[n,3,2] = inv_logit(beta_3 * x[n]);
         phi[n,3,3] = 1 -  phi[n,3,2]; 
  
         phi[n,2] = softmax(beta_2 * x[n]);

      for (i in 1:K){
        for (j in 1:K){
          Mmatrix[n,j,i] = phi[n,i,j];
        }
      }        
      Ps[n] = matrix_pow(Mmatrix[n], (period[n]-1)) * Ps_init;
  }

  for(n in 1:N){
      y[n] ~ categorical(Ps[n]);
  }
}



