require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./sim1_setup.R'); source("./sim1_functions.R")

tau.set = c(0.05,0.1,0.2)
bar.beta.set = c(0.0075,0.01,0.0125)

ss = matrix(c(84, 68, 62,
              52, 43, 41,
              36, 32, 31), nrow = 3, byrow = TRUE)

power = ss * 0

treatment.data = potential.effects(P, window.length)
treatment.data.weekend = potential.effects(P.weekend, window.length)

for(i in 1:length(bar.beta.set)) {
    for(j in 1:length(tau.set)) {
        tau = rep(tau.set[j],length(N))

        ## Treatment vector
        Z.t = Vectorize(cov.gen)((1:num.days) * T)
        d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
        daily.treat = -t(Z.t)%*%d

        num.persons = ss[i,j]
        num.iters = 1000

        initial.study = foreach(k=1:num.iters,
            .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
            estimation.simulation(num.persons, N, pi, tau, P, P.weekend,
                                  daily.treat, T, window.length, min.p, max.p,
                                  treatment.data, treatment.data.weekend)

        current.result = c(bar.beta.set[i], tau.set[j],ss[i,j], mean(initial.study))

        print(current.result)
        power[i,j] = mean(initial.study)

    }
}

save(results,file="sim1_result.RData")

stopCluster(cl)

