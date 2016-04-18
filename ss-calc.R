source("./functions.R"); source("./setup.R")

num.iters = 2000
x = 1
sim.output1 = foreach(i=1:num.iters, .combine = "cbind") %do% rho.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
x = 2
sim.output2 = foreach(i=1:num.iters, .combine = "cbind") %do% rho.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)

mean.rho.t1 = foreach(i=1:720, .combine = "c") %do% mean(sim.output1[i,sim.output1[i,]!=0])
mean.rho.t2 = foreach(i=1:720, .combine = "c") %do% mean(sim.output2[i,sim.output2[i,]!=0])

invvar.rho.t1 = foreach(i=1:720, .combine = "c") %do% mean(1/(sim.output1[i,sim.output1[i,]!=0]*(1-sim.output1[i,sim.output1[i,]!=0])))
invvar.rho.t2 = foreach(i=1:720, .combine = "c") %do% mean(1/(sim.output2[i,sim.output2[i,]!=0]*(1-sim.output2[i,sim.output2[i,]!=0])))

par(mfrow = c(2,1), mar = c(2,4,1,1)+0.1)
plot(1:720/60,mean.rho.t1, xlab = "", ylab = "Rho | S")
plot(1:720/60,mean.rho.t2, xlab = "Time (in hours)", ylab = "Rho | NS")

# smooth.rho1 = lowess(invvar.rho.t1,f=1/5)$y
# smooth.rho2 = lowess(invvar.rho.t2,f=1/5)$y

par(mfrow = c(2,1), mar = c(2,4,1,1)+0.1)
plot(1:720/60,invvar.rho.t1, xlab = "", ylab = "f | S")
plot(1:720/60,invvar.rho.t2, xlab = "Time (in hours)", ylab = "f | NS")

Z.t = Vectorize(cov.gen)(1:(days*T))

term1 = (rep(invvar.rho.t1,days))#*tau[1]*pi[1])
middle.term1 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term1[i]*outer(Z.t[,i],Z.t[,i]))

term2 = (rep(invvar.rho.t2,days))#*tau[2]*pi[2])
middle.term2 = foreach(i=1:ncol(Z.t), .combine = "+") %do% (term2[i]*outer(Z.t[,i],Z.t[,i]))

outer.term = foreach(j=1:ncol(Z.t), .combine = "+") %do% outer(Z.t[,j],Z.t[,j])

Sigma.1 = solve(outer.term,middle.term1)%*%solve(outer.term)/(sum(pi[1]*tau[1]))
Sigma.2 = solve(outer.term,middle.term2)%*%solve(outer.term)/(sum(pi[2]*tau[2]))

### sigma^2_x calculation

# num.iters = 500
# x = 1
# sigma.output1 = foreach(i=1:num.iters, .combine = "cbind") %do% var.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
# x = 2
# sigma.output2 = foreach(i=1:num.iters, .combine = "cbind") %do% var.function(x,N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)

sigmasq.1 = 0.004321966 #mean(sigma.output1,na.rm = TRUE)
sigmasq.2 = 0.003752347 #mean(sigma.output2,na.rm = TRUE)

b1 = d # Unstandardized effect sizes
b2 = d # Unstandardized effect sizes

samp.size.const = b1%*%solve(sigma.1* Sigma.1, b1) + b2%*%solve(sigma.2 * Sigma.2, b2)

sample.size(samp.size.const,p = 6,q = 6)
c(std.beta*0.004,tau)
## I only want to test for "Stress"
sample.size(b1%*%solve(sigma.1*Sigma.1, b1),p = 3,q = 3)
sample.size(b2%*%solve(Sigma.2, b1),p = 3,q = 3)


# ### Calculating the standardized effects
# phi1 = sigma.1/(prod(pi)/60) # approximately for now the variance scaling
# phi2 = sigma.2/(prod(pi)/60) # approximately for now the variance scaling
# gamma = 0.10 # percent decrease that we want to see
# pi.new = c(pi[1]*(1-gamma),1-pi[1]*(1-gamma))
# std.effect1 = pi[1]*gamma/sqrt(phi1*prod(pi.new)/60)
# std.effect2 = pi[1]*gamma/sqrt(phi2*prod(pi.new)/60)
#
# t = seq(1,days*T)
#
# Z.t = Vectorize(cov.gen)(t)
#
# d1 = find.d(bar.d=std.effect1,init.d=0,max.d=7,Z.t)
# d2 = find.d(bar.d=std.effect2,init.d=0,max.d=7,Z.t)
#
