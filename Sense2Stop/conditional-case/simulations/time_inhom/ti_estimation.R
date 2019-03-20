##First read in the arguments listed at the command line
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

all_numpersons = c(50,67,127)
all_barbeta = c(0.03, 0.025, 0.02)

source("./ti_setup.R")
source("./ti_functions.R")

num.persons = all_numpersons[all_barbeta == barbeta]

print(barbeta)
print(num.persons)

power = rep(0,2)


library(doParallel)
library(doRNG)

cl <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
if (is.na(cl)) { 
  cl <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
} 


registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()


for (j in 1:2) {
  ## Weekend transition matrix
  P.wkend = Pwkend.list[[j]]
  pi.wkend = piwkend.list[[j]]
  
  ## Treatment vector
  Z.t = Vectorize(cov.gen)((1:num.days) * T)
  d = find.d(barbeta,init.d,max.d,Z.t,num.days)
  daily.treat = -t(Z.t)%*%d
  
  P.treat.list = list()
  init.values.wkend = init.inputs.wkend[[j]]
  
  for(day in 1:num.days) {
    effect = rep(daily.treat[day],2)
    if((day == 6) | (day == 7)) {
      temp = optim(init.inputs.wkday,
                   effect.gap(P.wkend,
                              window.length, effect))
      P.treat.list[[day]] = calculateP(temp$par)
    } else {
      temp = optim(init.values.wkend,
                   effect.gap(P, window.length, effect))
      P.treat.list[[day]] = calculateP(temp$par)
    }
  }
  
  All.studies = foreach(k=1:1000, .combine = c, .options.RNG=231310) %dorng%
    estimation.simulation(num.persons, N, pi, P, pi.wkend, P.wkend,
                          P.treat.list, T, window.length, min.p, max.p)

  power[j] = mean(All.studies)
  
  print(c(barbeta, j, mean(All.studies)))
}

saveRDS(power, file = paste("power_",barbeta,".rds", sep = ""))

