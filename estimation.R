require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source("./functions.R")
bar.d = c(0.005,0.01,0.015,0.02); tau = c(0.1)
ss = c(123,37,24,19)+3
pc = matrix(nrow = length(bar.d), ncol = length(tau))

for(k in 1:length(tau)) {
  for(k.prime in 1:length(bar.d)) {
    
    current.tau = rep(tau[k],length(N))
    ### Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(bar.d[k.prime],init.d,max.d, Z.t,num.days)
    daily.treat = -t(Z.t)%*%d

    num.persons = ss[k.prime]
    
    print(c(bar.d[k.prime],current.tau[1],num.persons))
    
    num.iters = 1000
    
    initial.study = foreach(i=1:num.iters, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar% 
      estimation.simulation(num.persons, N, pi, current.tau, P.0, daily.treat, T, window.length, min.p, max.p)
    
    pc[k.prime,k] = mean(initial.study)
    
  }
}

pc

stopCluster(cl)
#save(initial.study,file="init_study.RData")
#write.table(ss,file="sample_size.txt")
write.table(pc,file="sample_size.txt")
#save(pc,file="pc.RData")

