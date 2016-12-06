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
bar.beta.set = c(0.0075,0.01,0.0125)

ss = power = matrix(0, nrow = length(bar.beta.set), ncol = length(tau.set))

treatment.data = potential.effects(P, window.length)

results = rep(0,0)

for(i in 1:length(bar.beta.set)) {
  for(j in 1:length(tau.set)) {
    tau = rep(tau.set[j],length(N))

    ### Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
    daily.treat = -t(Z.t)%*%d
    
    # Calculate Sample Size
    num.iters.ss = 500
    Sigma.params = ss.parameters(num.iters.ss, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)
    Q = Sigma.params[1:6,]; W = Sigma.params[7:12,]
    
    Sigma = solve(Q,W)%*%solve(Q)
    
    b1 =  d # Unstandardized effect sizes
    b2 = d # Unstandardized effect sizes
    beta = c(b1,b2)
    
    samp.size.const = beta%*%solve(Sigma, beta)
    
    initial.N = sample.size(samp.size.const,p = 6,q = 6)
    
    final.output = binary.search(initial.N, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
      
    ss[i,j] = final.output[1]
    power[i,j] = final.output[2]
  }
}

save(ss,file="ss.RData")
save(power,file="power.RData")

stopCluster(cl)

