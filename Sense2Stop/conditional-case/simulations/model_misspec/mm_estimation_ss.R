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
init.ss = c(118, 66, 41)
which.worst = c(15,7,7)

new.ss = new.power = matrix(nrow = length(bar.beta.set), ncol = 16)

for(i in 1:length(bar.beta.set)) {
    ## Baseline transition matrix
    P = P.list[[which.worst[i]]]
    pi = pi.list[[which.worst[i]]]

    ## Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
    daily.treat = -t(Z.t)%*%d

    P.treat.list = list()

    for(day in 1:num.days) {
        effect = rep(daily.treat[day],2)
        temp = optim(init.inputs,effect.gap(P, window.length, effect))
        P.treat.list[[day]] = calculateP(temp$par)
    }

    ## Initial sample size from prior simulations
    initial.N = init.ss[i]

    print(paste("Initial N =",initial.N))

    High.N.current = High.N.old = round(initial.N*1.15,0)
    Mid.N.current = Mid.N.old = initial.N
    Low.N.current = Low.N.old = round(initial.N*0.85,0)

    which.run = rep(TRUE, 3); power.old = power.current = rep(0,3)
    total.evals = 0
    num.old = c(Low.N.old,Mid.N.old,High.N.old)

    binary.iters = 100

    for(iter in 1:binary.iters) {
        print(paste("Iteration number",iter))
        High.N.old = High.N.current
        Mid.N.old = Mid.N.current
        Low.N.old = Low.N.current

        if(which.run[1] == TRUE) {
            set.seed("231310")
            Low.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
            estimation.simulation(Low.N.old, N, pi, P, P.treat.list, T, window.length, min.p, max.p)
            power.L.old = mean(Low.study)
        } else {power.L.old = power.L.current}

        set.seed("231310")
        Mid.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
        estimation.simulation(Mid.N.old, N, pi, P, P.treat.list, T, window.length, min.p, max.p)
        power.M.old = mean(Mid.study)

        if(which.run[3] == TRUE) {
            set.seed("231310")
            High.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dorng%
            estimation.simulation(High.N.old, N, pi, P, P.treat.list, T, window.length, min.p, max.p)
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

    new.ss[i] = final.output[1]
    new.power[i] = final.output[2]
}

print(max.value)

save(new.power,file="newpower.RData")
save(new.ss,file="newss.RData")

stopCluster(cl)

