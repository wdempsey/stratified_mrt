source("./functions.R")
# Study assumptions
num.persons = 17 # Number of persons in the study

T = 600 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

# Construct Transition Matrix given Inputs 1 and 2
bar.W = c(0.067, 0.519)
bar.Z = c(10.9,12.0)

tilde.Z = c((bar.Z[1]-3)/2, 0,(bar.Z[1]-3)/2,
(bar.Z[2]-3)/2, 0,(bar.Z[2]-3)/2)

P  = matrix(0, nrow = 6, ncol = 6)

diag(P) = tilde.Z/(1+tilde.Z)
P[2,3] = P[5,6] = 1.0
P[1,2] = 1-P[1,1]; P[4,5] = 1-P[4,4]
P[3,1] = (1-bar.W[1]) * (1-P[3,3]); P[3,4] = bar.W[1] * (1-P[3,3])
P[6,1] = (1-bar.W[2]) * (1-P[6,6]); P[6,4] = bar.W[2] * (1-P[6,6])

eig.P = eigen(P)

pi = (eig.P$vectors%*%diag(c(1,rep(0,5)))%*%solve(eig.P$vectors))[1,]  # Stationary distribution for P

# tau = rep(0.1,2) # Expected availability in each group ("Stressed", "Not Stressed")

# Randomization probability inputs
N =  c(0,1.61,0,0,2.21,0) # Avg. number of actions per day in each group ("Stressed", "Not Stressed")

lambda = 0.3 # Smoothing parameter in the randomization formula

min.p = 0.001 # Minimum randomization probability at each time point
max.p = 0.999 # Maximum randomization probability at each time point

# Treatment assumptions
init.d = 0 # Initial treatment effect
max.d = num.days/2 # Day of maximum treatment effect
bar.d = 0.01 # Avg treatment effect

## Initial inputs for finding the closest P.treat to give you the
## correct treatment effect.
init.inputs = c(P[1,1],P[3,3],P[3,4]/(1-P[3,3]), P[4,4], P[6,6], P[6,4]/(1-P[6,6]))

## Test that the choice of N leads to approximately 1.5 interventions per day 
## given stressed and not stressed under the null model.
# set.seed("81740")
# num.iters = 1000
# res.nonstress = res.stress = vector(length = num.iters)
# for (i in 1:num.iters) {
#    test = daily.sim(N,pi,P,P,T,window.length,min.p,max.p)
#    res.nonstress[i] = sum(test$A[test$X == 2])
#    res.stress[i] = sum(test$A[test$X == 5])
# }
# 
# mean(res.nonstress)
# mean(res.stress)
