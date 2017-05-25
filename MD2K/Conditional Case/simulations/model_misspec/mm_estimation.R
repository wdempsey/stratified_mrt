library(Rmpi)
library(parallel)
library(snow)
library(doParallel)

# reads the list of nodes which have been allocated
# by the cluster queue manager
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile, skip=1, header=FALSE)

ncpu <- mpi.universe.size() - 1

# builds a socket cluster using these nodes
cl <- makeCluster(ncpu, type='MPI')

print(ncpu)

registerDoParallel(cl)

library(doRNG)

source('./mm_setup.R'); source("./mm_functions.R")

bar.beta.set = c(0.02,0.025,0.03)
ss = c(118, 66, 41)

power = matrix(nrow = length(bar.beta.set), ncol = 16)
max.value = 0.0

for(i in 1:length(bar.beta.set)) {
    for (j in 1:16) {
        ## Number of persons given bar.beta.set[i]
        num.persons = ss[i]
        ## Baseline transition matrix
        P = P.list[[j]]
        pi = pi.list[[j]]

        ## Treatment vector
        Z.t = Vectorize(cov.gen)((1:num.days) * T)
        d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
        daily.treat = -t(Z.t)%*%d

        P.treat.list = list()

        for(day in 1:num.days) {
            effect = rep(daily.treat[day],2)
            temp = optim(init.inputs,effect.gap(P, window.length, effect))
            max.value = max(max.value, temp$value)
            P.treat.list[[day]] = calculateP(temp$par)
        }

        set.seed("231310")
        All.studies = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
            estimation.simulation(num.persons, N, pi, P, P.treat.list, T, window.length, min.p, max.p)

        power[i,j] = mean(All.studies)

        print(c(bar.beta.set[i], j, mean(All.studies)))
    }
}

print(max.value)

save(power,file="power.RData")

stopCluster(cl)

