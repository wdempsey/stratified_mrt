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

### Parallelize the optimizations!!!! S
### So all 10*3 = 30 cores needed for 6 - 8 hours!! 

bar.beta.set = c(0.02,0.0250,0.030)

Delta = window.length
output = p_all.k(Delta, theta.0)
baseline.prox = proximal.outcome(output, theta.0)
init.theta = unlist(theta.0)

results = optimal.treatment.barbetaset(baseline.prox, Delta, bar.beta.set, init.theta)

list.results = list()

for (i in 1:length(bar.beta.set)) {
  start = (i-1)*num.days + 1
  end = (i-1)*num.days + num.days
  list.results[[i]] = results[,start:end]
}


save(list.results,file="results.RData")
# save(power,file="power.RData")

stopCluster(cl)

