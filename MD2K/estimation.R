require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source("./functions.R")

#tau.set = c(0.05,0.1,0.2)
tau.set = c(0.1)
bar.beta.set = c(0.0075)
#bar.beta.set = c(0.0075,0.01,0.0125)
# ss = matrix(c(216,188,180,58,53,50,33,29,27,23,22,21), nrow = 4, byrow = TRUE)

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
    num.iters.ss = 1000
    Sigma.params = ss.parameters(num.iters.ss, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)
    Q = Sigma.params[1:6,]; W = Sigma.params[7:12,]
    
    Sigma = solve(Q,W)%*%solve(Q)
    
    b1 =  d # Unstandardized effect sizes
    b2 = d # Unstandardized effect sizes
    beta = c(b1,b2)
    
    samp.size.const = beta%*%solve(Sigma, beta)
    
    num.persons = sample.size(samp.size.const,p = 6,q = 6)
    
    poss.persons = num.persons+seq(-5,5,1)

    num.iters = 1

    for(peeps in poss.persons) {

        initial.study = foreach(k=1:num.iters, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
            estimation.simulation(peeps, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)

        current.result = c(bar.beta.set[i], tau.set[j],est_bar.d,peeps,mean(initial.study))

        print(current.result)
        result = c(result,current.result)
    }
  }
}

save(result,file="result.RData")

stopCluster(cl)

