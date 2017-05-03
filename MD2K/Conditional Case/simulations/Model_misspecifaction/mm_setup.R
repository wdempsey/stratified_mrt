# Study assumptions

num.persons = 17 # Number of persons in the study

T = 600 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

P = matrix(c(0.667,0.053,0.289,0.325), nrow = 2,byrow = TRUE)
#P = matrix(c(0.85,0.15,0.47,0.53), nrow = 2,byrow = TRUE)

P[1,] = P[1,]/sum(P[1,]); P[2,] = P[2,]/sum(P[2,])

eig.P = eigen(P)

pi = (eig.P$vectors%*%diag(c(1,0))%*%solve(eig.P$vectors))[1,]  # Stationary distribution for P

# tau = rep(0.1,2) # Expected availability in each group ("Stressed", "Not Stressed")

# Randomization probability inputs
N =  c(1.5,1.5) # Avg. number of actions per day in each group ("Stressed", "Not Stressed")

lambda = 0.3 # Smoothing parameter in the randomization formula

min.p = 0.001 # Minimum randomization probability at each time point
max.p = 0.999 # Maximum randomization probability at each time point

# Treatment assumptions
init.d = 0 # Initial treatment effect
max.d = num.days/2 # Day of maximum treatment effect
# bar.d = 0.01 # Avg treatment effect

#### Construct daily treatment vector
# source("./sim4_functions.R")
# Z.t = Vectorize(cov.gen)((1:num.days) * T)
# d = find.d(bar.d,init.d,max.d, Z.t,num.days)
# daily.treat = -t(Z.t)%*%d
# 
# 
