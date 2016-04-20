require(doParallel)

# reads the list of nodes which have been allocated 
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source('./ss-calc.R')

### Estimation procedure given the dataset
num.iters = 1000

initial.study = foreach(i=1:num.iters, .combine = c) %dopar% estimation.simulation(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

mean(initial.study)

stopCluster(cl)
save(initial.study,file="init_study.RData")


