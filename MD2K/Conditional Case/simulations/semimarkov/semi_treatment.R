## This file computes the treatment
## 
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

print(ncpu)

registerDoParallel(cl)

library(doRNG)

source('./semi_setup.R'); source("./semi_functions.R")

# bar.beta.set = c(0.02,0.0250,0.030)

bar.beta.set = 0.0250

system.time(temp <- p_all.k(window.length, theta.0))
system.time(temp2 <- parallel.p_all.k(window.length, theta.0))

# for(i in 1:length(bar.beta.set)) {
#   ## Treatment vector
#   Z.t = Vectorize(cov.gen)((1:num.days) * T)
#   d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
#   daily.treat = -t(Z.t)%*%d
#   
#   day = 5
#   
#   alt.beta = rep(daily.treat[5], 2)
#   
#   
#   
# }

# save(ss,file="ss.RData")
# save(power,file="power.RData")

stopCluster(cl)

