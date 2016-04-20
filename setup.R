# Study assumptions

num.persons = 18 # Number of persons in the study

T = 720 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

P.treat = P.0 = P = matrix(c(0.6,0.4,0.04,0.96), byrow = TRUE, nrow = 2, ncol = 2) # Initialized values of transition matrix

eig.P = eigen(P)

pi = (eig.P$vectors%*%diag(c(1,0))%*%solve(eig.P$vectors))[1,]  # Stationary distribution for P

tau = rep(0.1,2)#c(0.05,0.05) # Expected availability in each group ("Stressed", "Not Stressed")

# Randomization probability inputs
N =  c(2,2) # Avg. number of actions per day in each group ("Stressed", "Not Stressed")

lambda = 0.3 # Smoothing parameter in the randomization formula

min.p = 0.001 # Minimum randomization probability at each time point
max.p = 0.999 # Maximum randomization probability at each time point

# Treatment assumptions
#percent.Delta = 0.2 # % decrease in proximal outcome we wish to detect
init.d = 0 # Initial treatment effect
max.d = num.days/2 # Day of maximum treatment effect
bar.d = 0.02 # Avg treatment effect

source("./functions.R")
### Treatment vector
Z.t = Vectorize(cov.gen)((1:num.days) * T)
d = find.d(bar.d,init.d,max.d, Z.t,num.days)
daily.treat = t(Z.t)%*%d*1.33 # Times 1.3 to switch from the treatment scale to the effect of stationary distribution scale
