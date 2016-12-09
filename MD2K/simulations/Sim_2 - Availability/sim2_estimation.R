require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source("./sim3_setup.R"); source("./sim3_functions.R")

tau.set = c(0.05,0.1,0.2)
bar.beta.set = c(0.005,0.01,0.015,0.02)
# ss = matrix(c(216,188,180,58,53,50,33,29,27,23,22,21), nrow = 4, byrow = TRUE)

result1 = result2 = result3 = rep(0,0)

for(i in 1:length(bar.beta.set)) {
    for(j in 1:length(tau.set)) {
        mean.tau = tau.set[j]
        weekend.tau = mean.tau*1.2
        weekend = which(1:num.days%%7 == 0 | 1:num.days%%7 == 6)
        weekday.tau = (mean.tau-weekend.tau*length(weekend)/num.days)*(1-length(weekend)/num.days)^(-1)
        tau1.daily =  rep(weekday.tau,
            num.days) # Week vs Weekend availability in each group ("Stressed", "Not Stressed")
        tau1.daily[weekend] = weekend.tau

        init.tau = mean.tau*0.8; end.tau = mean.tau*1.2
        tau2.daily = seq(init.tau,end.tau, length.out = num.days) # Increasing availability over time (linear)
        tau3.daily = seq(end.tau,init.tau, length.out = num.days) # Decreasing availability over time (linear)

        ## Treatment vector
        Z.t = Vectorize(cov.gen)((1:num.days) * T)
        d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
        daily.treat = -t(Z.t)%*%d
        num.persons = ss[i,j]

        num.iters = 1000

        study1 = foreach(k=1:num.iters,
            .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
            estimation.simulation(num.persons, N, pi, tau1.daily, P, daily.treat,
                                  T, window.length, min.p, max.p)

        current.result1 = c(bar.beta.set[i], tau.set[j],mean(study1))

        print(current.result1)
        result1 = c(result1,current.result1)


        study2 = foreach(k=1:num.iters,
            .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
            estimation.simulation(num.persons, N, pi, tau2.daily, P, daily.treat,
                                  T, window.length, min.p, max.p)

        current.result2 = c(bar.beta.set[i], tau.set[j],mean(study2))

        print(current.result2)
        result2 = c(result2,current.result2)

        study3 = foreach(k=1:num.iters,
            .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
            estimation.simulation(num.persons, N, pi, tau3.daily, P, daily.treat,
                                  T, window.length, min.p, max.p)

        current.result3 = c(bar.beta.set[i], tau.set[j],mean(study3))

        print(current.result3)
        result3 = c(result3,current.result3)
    }
  }
}

save(result1,file="result1.RData")
save(result2,file="result2.RData")
save(result3,file="result3.RData")

stopCluster(cl)

