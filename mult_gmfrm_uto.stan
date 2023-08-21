data{
  int <lower=0> J; 
  int <lower=0> I; 
  int <lower=0> R; 
  int <lower=2> K; 
  int <lower=0> D; 
  int <lower=0> N; 
  int <lower=1, upper=J> ExamineeID [N]; 
  int <lower=1, upper=I> ItemID [N]; 
  int <lower=1, upper=R> RaterID [N]; 
  int <lower=1, upper=K> X [N]; 
}
transformed data{
  vector[K] c = cumulative_sum(rep_vector(1, K)) - 1;
}
parameters {
  vector[D] theta[J];
  real<lower=0> alpha_r [R-1];
  vector<lower=0>[D] alpha_i [I];
  vector[R-1] beta_r;
  vector[I] beta_i;
  vector[K-2] beta_ik [I];
}
transformed parameters{
  real<lower=0> trans_alpha_r[R];
  vector[R] trans_beta_r;
  vector[K-1] category_est[I];
  vector[K] category_prm[I];
  trans_alpha_r[1] = 1.0 / prod(alpha_r);
  trans_beta_r[1] = -1*sum(beta_r);
  trans_alpha_r[2:R] = alpha_r;  
  trans_beta_r[2:R] = beta_r;
  for(z in 1:I){
    category_est[z, 1:(K-2)] = beta_ik [z];
    category_est[z, K-1] = -1*sum(beta_ik [z]);  
    category_prm[z] = cumulative_sum(append_row(0, category_est[z]));
  }
}
model{
  for (d in 1:D){
    theta[,d] ~ normal(0, 1);
    alpha_i[,d] ~ lognormal(0, 1.0);
  }
  trans_alpha_r ~ lognormal(0, 1.0);
  trans_beta_r ~ normal(0, 1);
  beta_i ~ normal(0, 1);
  for (z in 1:I) category_est [z,] ~ normal(0, 1);
  for (n in 1:N){
    X[n] ~ categorical_logit(1.7 * trans_alpha_r[RaterID[n]] * (c*(dot_product(alpha_i[ItemID[n]], theta[ExamineeID[n]])-beta_i[ItemID[n]]-trans_beta_r[RaterID[n]])-category_prm[ItemID[n]]));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], 1.7 * trans_alpha_r[RaterID[n]] * (c*(dot_product(alpha_i[ItemID[n]], theta[ExamineeID[n]])-beta_i[ItemID[n]]-trans_beta_r[RaterID[n]])-category_prm[ItemID[n]]));
  }
}
