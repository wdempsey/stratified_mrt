##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

print(args)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  barbeta = 0.02
  day = 1
}else{
  barbeta = as.numeric(args[[1]])
  day = as.numeric(args[[2]])
}

print(barbeta)
print(day)

library(Rmpi)
library(parallel)
library(snow)
library(doParallel)

source('./semi_setup.R'); source("./semi_functions.R")

Delta = window.length
output = p_all.k(Delta, theta.0)
baseline.prox = proximal.outcome(output, theta.0)
init.theta = unlist(theta.0)

Z.t = Vectorize(cov.gen)((1:num.days) * T)
d = find.d(barbeta,init.d,max.d,Z.t,num.days)
daily.treat = -t(Z.t)%*%d

print(barbeta)
print(day)

results = optimal.treatment.day(baseline.prox, Delta, daily.treat, day, init.theta)

print(results)
