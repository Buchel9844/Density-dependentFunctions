// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;  // Number of plots/obs
  int<lower = 1> S;  // Number of species
  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
   int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
}

parameters{
  real<lower=0> lambdas[1];
  vector<lower=0>[S] m; //slope of the inflation point - specific distribution
  vector[S] c; //slope of the inflation point - specific distribution

  vector[S] alpha_function_tilde;
    
  real<lower=0> disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{
 // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  
 // loop parameters
  matrix[N,S] alpha_function_eij;
  vector[N] lambda_ei;
  
  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1];

    for(s in 1:S){
  
      alpha_function_eij[i,s]= alpha_function_tilde[s]*SpMatrix[i,s];

    }
  
    
 F_hat[i] = lambda_ei[i] / 1 + sum(alpha_function_eij[i,]);

  }

}

model{
  // set regular priors
  
  alpha_function_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);

 for(i in 1:N){
  Fecundity[i] ~ poisson(F_hat[i]); 
   }

}
generated quantities{
  vector[N] F_sim;
    if(run_estimation==1){
 for(i in 1:N){
    if(Fecundity[i] <= 0) break ;
    F_sim[i] = poisson_rng(Fecundity[i]);
              }
    }
}



