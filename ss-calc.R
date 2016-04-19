source("./functions.R"); source("./setup.R")

num.iters = 2000
x = 1
sim.output1 = foreach(i=1:num.iters, .combine = "cbind") %dopar% rho.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
x = 2
sim.output2 = foreach(i=1:num.iters, .combine = "cbind") %dopar% rho.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)

mean.rho.t1 = foreach(i=1:720, .combine = "c") %dopar% mean(sim.output1[i,sim.output1[i,]!=0])
mean.rho.t2 = foreach(i=1:720, .combine = "c") %dopar% mean(sim.output2[i,sim.output2[i,]!=0])

invvar.rho.t1 = foreach(i=1:720, .combine = "c") %dopar% mean(1/(sim.output1[i,sim.output1[i,]!=0]*(1-sim.output1[i,sim.output1[i,]!=0])))
invvar.rho.t2 = foreach(i=1:720, .combine = "c") %dopar% mean(1/(sim.output2[i,sim.output2[i,]!=0]*(1-sim.output2[i,sim.output2[i,]!=0])))

par(mfrow = c(2,1), mar = c(2,4,1,1)+0.1)
plot(1:720/60,mean.rho.t1, xlab = "", ylab = "Rho | S")
plot(1:720/60,mean.rho.t2, xlab = "Time (in hours)", ylab = "Rho | NS")

par(mfrow = c(2,1), mar = c(2,4,1,1)+0.1)
plot(1:720/60,invvar.rho.t1, xlab = "", ylab = "f | S")
plot(1:720/60,invvar.rho.t2, xlab = "Time (in hours)", ylab = "f | NS")

## Fix NaNs by taking average of the values (for now)
invvar.rho.t1[is.nan(invvar.rho.t1)] = mean(invvar.rho.t1,na.rm=TRUE)
invvar.rho.t2[is.nan(invvar.rho.t2)] = mean(invvar.rho.t2,na.rm=TRUE)

Z.t = Vectorize(cov.gen)(1:(num.days*T))
d = find.d(bar.d, init.d, max.d, Z.t, num.days)

term1 = (rep(invvar.rho.t1,num.days))#*tau[1]*pi[1])
middle.term1 = foreach(i=1:ncol(Z.t), .combine = "+") %dopar% (term1[i]*outer(Z.t[,i],Z.t[,i]))

term2 = (rep(invvar.rho.t2,num.days))#*tau[2]*pi[2])
middle.term2 = foreach(i=1:ncol(Z.t), .combine = "+") %dopar% (term2[i]*outer(Z.t[,i],Z.t[,i]))

outer.term = foreach(j=1:ncol(Z.t), .combine = "+") %dopar% outer(Z.t[,j],Z.t[,j])

Sigma.1 = solve(outer.term,middle.term1)%*%solve(outer.term)/(sum(pi[1]*tau[1]))
Sigma.2 = solve(outer.term,middle.term2)%*%solve(outer.term)/(sum(pi[2]*tau[2]))

### sigma^2_x calculation

# num.iters = 500
# x = 1
# sigma.output1 = foreach(i=1:num.iters, .combine = "cbind") %dopar% var.function(x,N, pi, tau, P.0, P.treat, T, windoparw.length, min.p, max.p)
# x = 2
# sigma.output2 = foreach(i=1:num.iters, .combine = "cbind") %dopar% var.function(x,N, pi, tau, P.0, P.treat, T, windoparw.length, min.p, max.p)

sigmasq.1 = 0.004321966 #mean(sigma.output1,na.rm = TRUE)
sigmasq.2 = 0.003752347 #mean(sigma.output2,na.rm = TRUE)

b1 = d # Unstandardized effect sizes
b2 = d # Unstandardized effect sizes

samp.size.const = b1%*%solve(sigmasq.1* Sigma.1, b1) + b2%*%solve(sigmasq.2 * Sigma.2, b2)

ss.answer = sample.size(samp.size.const,p = 6,q = 6)

ss.answer
