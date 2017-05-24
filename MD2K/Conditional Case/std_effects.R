source('./setup.R'); source("./functions.R")
library(doRNG)

bar.beta.set = c(0.02,0.025,0.03)

ss = power = vector(length = length(bar.beta.set))

for(i in 1:length(bar.beta.set)) {

    ### Treatment vector
    Z.t = Vectorize(cov.gen)((1:num.days) * T)
    d = find.d(bar.beta.set[i],init.d,max.d,Z.t,num.days)
    daily.treat = -t(Z.t)%*%d

    P.treat.list = list()

    for(day in 1:num.days) {
        effect = rep(daily.treat[day],2)
        temp = optim(init.inputs,effect.gap(P, window.length, effect))
        P.treat.list[[day]] = calculateP(temp$par)
    }

    ## Calculate Sample Size
    num.iters = 100
    test.sigma = barsigma.estimation(num.iters, N, pi, P, P.treat.list, T,
                                window.length, min.p, max.p)
    ## test.tildepr = tildepr.estimation(num.iters, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p)
    print(c(bar.beta.set[i], round(bar.beta.set[i]/sqrt(test.sigma[,1]/test.sigma[,2]),4)))
    ## print(c(bar.beta.set[i], tau.set[j], test.tildepr[,1]/test.tildepr[,2]))
  }
}
