source('./setup.R'); source("./functions.R")

tau.set = c(0.05,0.1,0.2)
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
    num.iters.ss = 250
    #test.sigma = barsigma.estimation(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)
    test.tildepr = tildepr.estimation(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)
    # print(c(bar.beta.set[i], tau.set[j], round(bar.beta.set[i]/sqrt(test.sigma[,1]/test.sigma[,2]),2)))
    print(c(bar.beta.set[i], tau.set[j], test.tildepr[,1]/test.tildepr[,2]))
  }
}