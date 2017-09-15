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
  numpersons = 17
}else{
  barbeta = as.numeric(args[[1]])
  numpersons = as.numeric(args[[2]])
}

print(barbeta)
print(numpersons)

all_treatmentthetas = read.csv("output/export_adjusted.csv", header = FALSE)

theta.treat.list = list()
for (day in 1:num.days) {
  if(day == 1) {
    theta.treat.list[[day]] = theta.0
  } else{
    row.id = which(all_treatmentthetas[,1] == -barbeta & all_treatmentthetas[,2] == day)
    unlisted.theta.day = as.numeric(all_treatmentthetas[row.id, 3:ncol(all_treatmentthetas)])
    theta.treat.list[[day]] = relist.thetas(unlisted.theta.day)
  }
}

library(Rmpi)
library(parallel)
library(snow)
library(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

ncpu <- mpi.universe.size() - 1

# builds a socket cluster using these nodes
cl <- makeCluster(ncpu, type='MPI')

source('./semi_setup.R'); source("./semi_functions.R")

Delta = window.length

set.seed("231310")
All.studies = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
  estimation.simulation(num.persons, N, pi, theta.0, theta.treat.list,
                        T, window.length, min.p, max.p, pi.SMC)

power = mean(All.studies)

print(c(barbeta, mean(All.studies)))

saveRDS(power, file = paste("power_",barbeta,".rds", sep = ""))

stopCluster(cl)

