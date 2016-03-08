### Randomization probabilities calculation
source("./functions.R")

## Initial conditions
T = 12 * 60/5 # number decision points per day, 12 hours possible every five minutes
days = 14 # number of days
N = c(2,2,0) # average number of interventions ("stressed", "not stressed", "unavailable")

## For now, I will explicitly define the transition matrix
P = matrix(c(0.403,0.271,0.327,0.029,0.704,0.268,0.063,0.362,0.575), nrow = 3, ncol = 3, byrow = TRUE)
M = eigen(P)$vectors
stationary = M%*%diag(c(1,0,0))%*%solve(M)
pi = stationary[1,]

### Availability is just sum of first two stationary components
min.p = min(N[1:2]/(T*pi[1:2]))/10  # Set to 1/10 of initial probability
max.p = 0.95
inv.var.est = inverse.expectation(1000,pi,P,T,N,min.p,max.p)
gamma = sum(pi[1:2])
#flat.est = c( quantile(inv.var.est[1,],probs=0.9),quantile(inv.var.est[2,],probs=0.9))
flat.est = rowMeans(inv.var.est)

### Initial sample size calculations
c(sample.size(gamma,inv.var.est, 0.10,T,days),
  sample.size(gamma,flat.est, 0.09,T,days),
  sample.size(gamma,flat.est, 0.08,T,days),
  sample.size(gamma,flat.est, 0.07,T,days),
  sample.size(gamma,flat.est, 0.06,T,days),
  sample.size(gamma,flat.est, 0.05,T,days)
)