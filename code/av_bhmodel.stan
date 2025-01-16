// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  vector[N] av;
  vector[N] mh;
  int<lower = 1>P;
  int Plot[N];
  real ag;
  real mg;
}

parameters{
  real epsilon[P];
  real<lower = 0> sigma;
  real<lower = 0> lambda;
  real alpha_av;
  real alpha_mh;

  //real<lower = 0>  alpha_av;
  //real<lower = 0>  alpha_mh;


}

model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma ~ gamma(0.001, 0.001);
  epsilon ~ gamma(sigma, sigma);
  alpha_av ~ normal(0, 10);
  alpha_mh ~ normal(0, 10);
  //lambda ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);



  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i]*ag / (1 + alpha_av*av[i]*ag + alpha_mh*mh[i]*mg);
    F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}





