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

source("./ti_setup.R")
source("./ti_functions.R")

bar.beta.set = c(0.02,0.025,0.03)
ss = c(114, 79, 54)

power = matrix(nrow = length(bar.beta.set), ncol = 16)

for(i in 1:length(bar.beta.set)) {
    for (j in 1:2) {
        ## Number of persons given bar.beta.set[i]
        num.persons = ss[i]

        ## Weekend transition matrix
        P.wkend = Pwkend.list[[j]]
        pi.wkend = piwkend.list[[j]]

        ## Treatment vector
        Z.t = Vectorize(cov.gen)((1:num.days) * T)
        d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
        daily.treat = -t(Z.t)%*%d

        P.treat.list = list()
        init.values.wkend = init.inputs.wkend[[j]]

        for(day in 1:num.days) {
            effect = rep(daily.treat[day],2)
            if((day == 6) | (day == 7)) {
                temp = optim(init.inputs.wkday,
                             effect.gap(P.wkend,
                                        window.length, effect))
                P.treat.list[[day]] = calculateP(temp$par)
            } else {
                temp = optim(init.values.wkend,
                             effect.gap(P, window.length, effect))
                P.treat.list[[day]] = calculateP(temp$par)
            }
        }

    set.seed("231310")
    All.studies = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
        estimation.simulation(num.persons, N, pi, P, pi.wkend, P.wkend,
                              P.treat.list, T, window.length, min.p, max.p)

    power[i,j] = mean(All.studies)

    print(c(bar.beta.set[i], j, mean(All.studies)))
  }
}

save(power,file="power.RData")

stopCluster(cl)
