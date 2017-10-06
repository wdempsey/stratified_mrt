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

all_numpersons = c(41, 66, 118)
all_barbeta = c(0.03, 0.025, 0.02)

num.persons = all_numpersons[all_barbeta == barbeta]

print(barbeta)
print(num.persons)

# library(Rmpi)
# library(parallel)
# library(snow)
library(foreach)
library(doRNG)
library(doParallel)

getDoParWorkers()

# Calculate the number of cores
no_cores <- detectCores() - 1

cl<-makeCluster(no_cores)
registerDoParallel(cl)

source('./semi_setup.R'); source("./semi_functions.R")

all_treatmentthetas = read.csv("output/export.csv", header = FALSE)

theta.treat.list = list()
for (temp.day in 1:num.days) {
  if(temp.day == 1) {
    theta.treat.list[[temp.day]] = theta.0
  } else{
    row.id = which(all_treatmentthetas[,1] == -barbeta & all_treatmentthetas[,2] == temp.day)
    unlisted.theta.day = as.numeric(all_treatmentthetas[row.id, 3:ncol(all_treatmentthetas)])
    theta.treat.list[[temp.day]] = relist.thetas(unlisted.theta.day)
  }
}

Delta = window.length


set.seed("231310")
All.studies = foreach(k=1:20, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
  estimation.simulation(num.persons, N, pi.simple, theta.0, theta.treat.list,
                        T, window.length, min.p, max.p, pi.SMC)

power = mean(All.studies)

print(c(barbeta, mean(All.studies)))

saveRDS(power, file = paste("power_",barbeta,".rds", sep = ""))

stopImplicitCluster()
