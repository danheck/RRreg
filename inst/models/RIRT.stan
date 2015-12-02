
# Stan RR IRT model

data{
  # hyperpriors
	real<lower=0> itemPrior;  # maxium bound for uniform prior on a and b
	
	# data and indices
	int<lower=1> N;  					# number of persons
	int<lower=1> K;  					# number of items
	int<lower=0, upper=1> Y[N, K];	    # chosen responses of partipants
	matrix<lower=0, upper=1>[K,2] P;
}

parameters{
  vector[N] theta;
  
  vector<lower=0, upper=itemPrior>[K] a;
  vector<lower=-itemPrior, upper=itemPrior>[K] b;
}


transformed parameters{
	// int<lower=0, upper=1> Y_true[N, K];	    # chosen responses of partipants
  matrix<lower=0, upper=1>[N,K] p;
  matrix<lower=0, upper=1>[N,K] pi;

  for(i in 1:N){
    for(k in 1:K){
      # IRT model: probability for true state
      pi[i,k] <- Phi_approx(a[k] * theta[i] - b[k]);
      
      # RR model: probability to respond 1=yes
      # (defined by misclassification matrix P)
      p[i,k] <- P[k,1]*(1-pi[i,k]) + P[k,2]*pi[i,k];
    }
  }

}

model{

  theta ~ normal(0, 1);
  
  for(i in 1:N){
    for(k in 1:K){
      # distribution of RR responses: model dependent
      Y[i,k] ~ bernoulli(p[i,k]);
//       Y_true[i,k] ~ bernoulli(pi[i,k]);
//       Y[i,k] ~ bernoulli(P[k,Y_true[k,i]+1]);
    }
  }
  
}