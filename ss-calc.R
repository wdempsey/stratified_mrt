source("./setup.R"); source("./functions.R")

num.iters = 200

treatment.data = potential.effects(P)

sim.output = foreach(i=1:num.iters, .combine = "rbind") %do% rho.function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)

sim.output = data.frame(sim.output)

# sim.output$middle = ((sim.output$A.t/sim.output$rho.t+(1-sim.output$A.t)/(1-sim.output$rho.t))*sim.output$psi.t^2)#^2
sim.output$middle = (1/(sim.output$rho.t*(1-sim.output$rho.t)))#*sim.output$psi.t#^2

# sim.aggregate = aggregate(middle~A.t+X.t+day,data = sim.output, mean, na.rm = TRUE)

### Variance Calculations
# var.output = foreach(i=1:num.iters, .combine = "rbind") %do% var.function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)

# sigmasq = aggregate(Y.t~A.t+X.t, data = var.output, mean, na.rm = TRUE)
sigmasq = 0.005#mean(var.output)

sim.aggregate$sigmasq = sim.aggregate$middle

### Aggregate over Actions
total.agg = aggregate(sigmasq~X.t+day,data=sim.aggregate, sum, na.rm = TRUE)

daily.sim.terms = total.agg$sigmasq

Z.t = Vectorize(cov.gen)(1:(num.days*T))
# d = find.d(bar.d, init.d, max.d, Z.t, num.days)

term1 = rep(daily.sim.terms[total.agg[,1]==1],each = T)
middle.term1 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term1[i]*outer(Z.t[,i],Z.t[,i]))

term2 = rep(daily.sim.terms[total.agg[,1]==2],each = T)
middle.term2 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term2[i]*outer(Z.t[,i],Z.t[,i]))

outer.term = foreach(j=1:ncol(Z.t), .combine = "+") %do% outer(Z.t[,j],Z.t[,j])

Sigma.1 = solve(outer.term,middle.term1)%*%solve(outer.term)/(sum(pi[1]*tau[1]))
Sigma.2 = solve(outer.term,middle.term2)%*%solve(outer.term)/(sum(pi[2]*tau[2]))

b1 =  d # Unstandardized effect sizes
b2 = d # Unstandardized effect sizes

samp.size.const = b1%*%solve(sigmasq*Sigma.1, b1) + b2%*%solve(sigmasq*Sigma.2, b2)

ss.answer = sample.size(samp.size.const,p = 6,q = 6)

paste("Sample size of ", ss.answer, " when percent change desired to detect is ", bar.d, ", with conditional expected availability both equal to ", tau[1],".", sep = "")
