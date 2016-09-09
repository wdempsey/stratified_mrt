require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source("./functions.R")

tau.set = c(0.05,0.1,0.2)
bar.beta.set = c(0.005,0.01,0.015,0.02)
# ss = matrix(c(216,188,180,58,53,50,33,29,27,23,22,21), nrow = 4, byrow = TRUE)

ss = pc = matrix(nrow = length(bar.beta.set), ncol = length(tau.set))

for(i in 1:length(bar.beta.set)) {
  for(j in 1:length(tau.set)) {
    tau = rep(tau.set[j],length(N))
    
    ### Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
    daily.treat = -t(Z.t)%*%d
    
    # Calculate the Sample Size --> Still and issue that the ss is highly variable
    num.iters.ss = 200
    Sigma.params = ss.parameters(num.iters.ss, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)
    Q = Sigma.params[1:6,]; W = Sigma.params[7:12,]
    bar.sigma.sq = 5.113 * 10^(-3)
    
    Sigma = solve(Q,W)%*%solve(Q)
    
    b1 =  d # Unstandardized effect sizes
    b2 = d # Unstandardized effect sizes
    beta = c(b1,b2)
    
    samp.size.const = beta%*%solve(bar.sigma.sq*Sigma, beta)
    
    num.persons = sample.size(samp.size.const,p = 6,q = 6)    
    
    ss[i,j] = num.persons
      
    print(c(bar.beta.set[i], tau.set[j],num.persons))
          
    num.iters = 1000
    
    initial.study = foreach(k=1:num.iters, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar% 
      estimation.simulation(num.persons, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)

    ss[i,j] = num.persons
    pc[i,j] = mean(initial.study)
    
    print(c(bar.beta.set[i], tau.set[j],num.persons,mean(initial.study)))
  }    
}

print(ss)
print(pc)

save(ss,file="sample_size.RData")
save(pc,file="power.RData")

stopCluster(cl)

