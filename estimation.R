source('./functions.R'); source('./setup.R')

### Estimation procedure given the dataset
Z.t = Vectorize(cov.gen)((1:num.days) * T)
d = find.d(bar.d,init.d,max.d, Z.t,num.days)
daily.treat = t(Z.t)%*%d

num.iters = 10

foreach(i=1:num.iters, .combine = c) %do% estimation.simulation(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)


system.time(estimation.simulation(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p))