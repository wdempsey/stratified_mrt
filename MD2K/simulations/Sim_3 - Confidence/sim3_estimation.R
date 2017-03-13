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

ss = matrix(c(84, 68, 62,
              52, 43, 41,
              36, 32, 31), nrow = 3, byrow = TRUE)

treatment.data = potential.effects(P, window.length)

C = c(0.01,0.1,1.0,10.0)

power = array(0, dim = c(length(bar.beta.set), length(tau.set), length(C)))

results = rep(0,0)

for(i in 1:length(bar.beta.set)) {
    for(j in 1:length(tau.set)) {
        for(k in 1:length(C)) {
            tau = rep(tau.set[j],length(N))

            ## Treatment vector
            Z.t = Vectorize(cov.gen)((1:num.days) * T)
            d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
            daily.treat = -t(Z.t)%*%d

            num.persons = ss[i,j]

            num.iters = 1000

            initial.study = foreach(k=1:num.iters,
                .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
                estimation.simulation(num.persons, N, pi, tau, C[k],
                                      P, daily.treat, T, window.length, 
                                      min.p, max.p, treatment.data)

            current.result = c(bar.beta.set[i], tau.set[j],C[k],mean(initial.study))
            
            power[i,j,k] = mean(initial.study)

            print(current.result)
            result = c(result,current.result)
        }
    }
}

save(result,file="result.RData")
save(power,file="power.RData")

stopCluster(cl)

