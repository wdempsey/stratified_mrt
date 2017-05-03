require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./setup.R'); source("./functions.R")

tau.set = c(1.00)
bar.beta.set = c(0.0075,0.01,0.0125)

ss = power = matrix(0, nrow = length(bar.beta.set), ncol = length(tau.set))

treatment.data = potential.effects(P, window.length)

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
    Q = Sigma.params[1:3,]; W = Sigma.params[4:6,]
    
    Sigma = solve(Q,W)%*%solve(Q)
    
    beta =  d # Unstandardized effect sizes
    
    samp.size.const = beta%*%solve(Sigma, beta)
    
    initial.N = sample.size(samp.size.const,p = 3,q = 3)
    
    print(paste("Initial N =",initial.N))
    
    High.N.current = High.N.old = round(initial.N*1.15,0)
    Mid.N.current = Mid.N.old = initial.N
    Low.N.current = Low.N.old = round(initial.N*0.85,0)
    
    which.run = rep(TRUE, 3); power.old = power.current = rep(0,3)
    total.evals = 0
    num.old = c(Low.N.old,Mid.N.old,High.N.old)
    
    binary.iters = 100
    
    for(iter in 1:binary.iters) {
      print(iter)
      High.N.old = High.N.current
      Mid.N.old = Mid.N.current
      Low.N.old = Low.N.current
      
      if(which.run[1] == TRUE) {
        Low.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
          estimation.simulation(Low.N.old, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
        power.L.old = mean(Low.study)
      } else {power.L.old = power.L.current}
      
      Mid.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
        estimation.simulation(Mid.N.old, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
      power.M.old = mean(Mid.study)
      
      if(which.run[3] == TRUE) {
        High.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
          estimation.simulation(High.N.old, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
        power.H.old = mean(High.study)
      } else {power.H.old = power.H.current}
      
      total.evals = total.evals + sum(which.run)
      
      num.old = c(Low.N.old,Mid.N.old,High.N.old)
      power.old = c(power.L.old,power.M.old,power.H.old)
      
      print(rbind(num.old,power.old))
      print(paste("Total Evaluations so far =", total.evals))
      
      if(power.L.old > power.M.old) {
        temp = mean(c(power.L.old,power.M.old))
        power.L.old = power.M.old = temp
      } 
      
      if(power.M.old > power.H.old) {
        temp = mean(c(power.H.old,power.M.old))
        power.H.old = power.M.old = temp
      } 
      
      if(power.L.old > power.H.old) {
        temp = mean(c(power.H.old,power.M.old,power.L.old))
        power.H.old = power.M.old = power.L.old = temp
      }
      
      power.fixed = c(power.L.old,power.M.old,power.H.old)
      
      if(High.N.old - Mid.N.old <= 1 & Mid.N.old - Low.N.old <= 1) {
        break
      }
      
      if(power.L.old > 0.80) {
        fit.power = lm(num.old~power.fixed)
        Mid.N.current = ceiling(fit.power$coefficients[1] + fit.power$coefficients[2]*0.80)
        Low.N.current = round(Mid.N.current*0.8,0)
        High.N.current = Low.N.old
        power.H.current = power.L.old
        which.run = c(TRUE, TRUE, FALSE)  
      } else if (power.H.old < 0.80) {
        fit.power = lm(num.old~power.fixed)
        Mid.N.current = ceiling(fit.power$coefficients[1] + fit.power$coefficients[2]*0.80)
        High.N.current = round(Mid.N.current*1.15,0)
        Low.N.current = High.N.old
        power.L.current = power.H.old
        which.run = c(FALSE, TRUE, TRUE)
      } else if (power.M.old < 0.80 ) {
        Mid.N.current = round(mean(c(Mid.N.old,High.N.old)),0)
        Low.N.current = Mid.N.old
        power.L.current = power.M.old
        which.run = c(FALSE, TRUE, TRUE)
      } else {
        Mid.N.current = round(mean(c(Mid.N.old,Low.N.old)),0)
        High.N.current = Mid.N.old
        power.H.current = power.M.old
        which.run = c(TRUE, TRUE, FALSE)
      }
      
    } 
    
    best = min(which(power.old>0.80))
    final.output = c(num.old[best], power.fixed[best])
    
    ss[i,j] = final.output[1]
    power[i,j] = final.output[2]
  }
}

save(ss,file="ss.RData")
save(power,file="power.RData")

stopCluster(cl)

