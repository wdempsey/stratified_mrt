source("./setup.R"); source("./functions.R")

num.persons = 100

ss.params = ss.parameters(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

Z.t = Vectorize(cov.gen)(1:(num.days*T))
# d = find.d(bar.d, init.d, max.d, Z.t, num.days)

term1 = rep(ss.params$exp.est[ss.params$exp.est[,1]==1,3],each = T)
middle.term1 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term1[i]*outer(Z.t[,i],Z.t[,i]))

term2 = rep(ss.params$exp.est[ss.params$exp.est[,1]==2,3],each = T)
middle.term2 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term2[i]*outer(Z.t[,i],Z.t[,i]))

outer.term = foreach(j=1:ncol(Z.t), .combine = "+") %do% outer(Z.t[,j],Z.t[,j])

Sigma.1 = solve(outer.term,middle.term1)%*%solve(outer.term)/(sum(pi[1]*tau[1]))
Sigma.2 = solve(outer.term,middle.term2)%*%solve(outer.term)/(sum(pi[2]*tau[2]))

b1 =  d # Unstandardized effect sizes
b2 = d # Unstandardized effect sizes

samp.size.const = b1%*%solve(ss.params$sigmasq*Sigma.1, b1) + b2%*%solve(ss.params$sigmasq*Sigma.2, b2)

ss.answer = sample.size(samp.size.const,p = 6,q = 6)

paste("Sample size of ", ss.answer, " when percent change desired to detect is ", bar.d, ", with conditional expected availability both equal to ", tau[1],".", sep = "")

#### Alternative including extra term

gamma = 0.8

alt.Sigma.1 = solve(outer.term,middle.term1*(1+tau[1]*sum(gamma^(1:60))))%*%solve(outer.term)/(sum(pi[1]*tau[1]))
alt.Sigma.2 = solve(outer.term,middle.term2*(1+tau[2]*sum(gamma^(1:60))))%*%solve(outer.term)/(sum(pi[2]*tau[2]))

alt.samp.size.const = b1%*%solve(ss.params$sigmasq*alt.Sigma.1, b1) + b2%*%solve(ss.params$sigmasq*alt.Sigma.2, b2)

alt.ss.answer = sample.size(alt.samp.size.const,p = 6,q = 6)

paste("Conservative sample size of ", alt.ss.answer, " when percent change desired to detect is ", bar.d, ", with conditional expected availability both equal to ", tau[1],".", sep = "")

