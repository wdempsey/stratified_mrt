require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source("./functions.R")
load("ss.RData")

pc = vector(length = nrow(ss.data))

for(k in 1:nrow(ss.data)) {
    current.tau = rep(ss.data[k,2],length(N))
    ### Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(ss.data[k,1],init.d,max.d, Z.t,num.days)
    daily.treat = -t(Z.t)%*%d

    num.persons = ss.data[k,3]
    
    num.iters = 1000
    
    initial.study = foreach(i=1:num.iters, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar% 
      estimation.simulation(num.persons, N, pi, current.tau, P.0, daily.treat, T, window.length, min.p, max.p)
    
    pc[k] = mean(initial.study)
    
    print(c(ss.data[k,],pc[k]))
    
}

ss.data = cbind(ss.data,pc)

stopCluster(cl)
save(ss.data,file="ss_with_power.RData")

