#First read in the arguments listed at the command line
args=(commandArgs(TRUE))

print(args)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  barbeta = 0.025
}else{
  barbeta = as.numeric(args[[1]])
}

barbeta

source('./mm_setup.R'); source("./mm_functions.R")

bar.beta.set = c(0.02,0.025,0.03)
ss.set = c(127, 67, 50)

## Number of persons given bar.beta.set[i]
num.persons = ss.set[bar.beta.set == barbeta]

power = 0.0
max.value = 0.0

library(doParallel)
library(doRNG)

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

for (j in 1:16) {
  ## Baseline transition matrix
  P = P.list[[j]]
  pi = pi.list[[j]]
  
  ## Treatment vector
  Z.t = Vectorize(cov.gen)((1:num.days) * T)
  d = find.d(barbeta,init.d,max.d,Z.t,num.days)
  daily.treat = -t(Z.t)%*%d
  
  P.treat.list = list()
  
  for(day in 1:num.days) {
    effect = rep(daily.treat[day],2)
    temp = optim(init.inputs,effect.gap(P, window.length, effect))
    max.value = max(max.value, temp$value)
    P.treat.list[[day]] = calculateP(temp$par)
  }
  
  All.studies = foreach(k=1:1000, .combine = c, .options.RNG=231310) %dorng%
    estimation.simulation(num.persons, N, pi, P, P.treat.list, T, window.length, min.p, max.p)
  
  power[j] = mean(All.studies)
  
  print(c(barbeta, j, mean(All.studies)))
}

print(max.value)

save(power,file="power.RData")

stopCluster(cl)

