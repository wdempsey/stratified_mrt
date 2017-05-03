require(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

# builds a socket cluster using these nodes
cl <- makeCluster(c(as.character(hostlist$V1)), type='SOCK')

registerDoParallel(cl)

source('./sim3_setup.R'); source("./sim3_functions.R")

tau.set = c(0.05,0.1,0.2)
bar.beta.set = c(0.0075,0.01,0.0125)

ss = matrix(c(71, 62, 58,
              44, 39, 37,
              30, 26, 25), nrow = 3, byrow = TRUE)

epsilon = c(0.02,0.002)

power = array(0, dim = c(length(bar.beta.set), length(tau.set), length(epsilon)))

result = rep(0,0)

P.epsilon = P

for(i in 1:length(bar.beta.set)) {
    for(j in 1:length(tau.set)) {
        for(k in 1:length(epsilon)) {
            tau = rep(tau.set[j],length(N))

            ## Treatment vector
            Z.t = Vectorize(cov.gen)((1:num.days) * T)
            d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
            daily.treat = -t(Z.t)%*%d

            num.persons = ss[i,j]

            num.iters = 1000
            
            P.epsilon[1,1] = P[1,1] - epsilon[k]
            P.epsilon[2,2] = P[2,2] - epsilon[k]
            P.epsilon[1,2] = P[1,2] + epsilon[k]
            P.epsilon[2,1] = P[2,1] + epsilon[k]
            
            treatment.data.epsilon = potential.effects(P.epsilon,window.length)

            study.epsilon = foreach(k=1:num.iters,
                .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
                estimation.simulation(num.persons, N, pi, tau,
                                      P.epsilon, daily.treat, T, window.length,
                                      min.p, max.p, treatment.data.epsilon)

            current.result = c(bar.beta.set[i], tau.set[j],epsilon[k],mean(study.epsilon))

            power[i,j,k] = mean(study.epsilon)

            print(current.result)
            result = c(result,current.result)
        }
    }
}

save(result,file="result.RData")
save(power,file="power.RData")

stopCluster(cl)

