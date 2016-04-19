source('./functions.R'); source('./setup.R')

### Estimation procedure given the dataset
Z.t = Vectorize(cov.gen)((1:num.days) * T)
d = find.d(bar.d,init.d,max.d, Z.t,num.days)
daily.treat = t(Z.t)%*%d

people = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

colnames(people) = c("person", "day","t","Y.t","A.t","X.t","rho.t","I.t")

Y.t.person = people[,4]
A.t.person = people[,5]
X.t.person = people[,6]
rho.t.person = people[,7]

B.t.person = t(Vectorize(cov.gen)(people[,2]*people[,3]))

rho = mean(mean.rho.t1)*pi[1]+mean(mean.rho.t2)*pi[2] # Need to think about choice of rho

Z.t.person = B.t.person*matrix(rep(people[,4]-rho,3), byrow=TRUE, ncol = 3)

cov.t.person = cbind(B.t.person,Z.t.person)

log.weights = A.t.person*(log(rho) - log(rho.t.person)) + (1-A.t.person)*(log(1-rho) - log(1-rho.t.person))

fit.people = lm(Y.t.person~cov.t.person:as.factor(X.t.person)-1,weights = exp(log.weights))
fit.people2 = lm(Y.t.person~B.t.person:as.factor(X.t.person)-1,weights = exp(log.weights))

Covariates = model.matrix(fit.people)

XWX = foreach(i=1:num.persons, .combine = "+") %do% extract.tXWX(Covariates,people,log.weights,i)

Middle = foreach(person=1:num.persons, .combine = "+") %do% M.function(Covariates, people, log.weights, person,XWX,fit.people)

Sigma = solve(XWX,Middle)%*%solve(XWX)

entries = c(4:6,10:12)
output = (fit.people$coefficients[entries]%*%solve(Sigma[entries,entries], fit.people$coefficients[entries]))/num.persons

alpha.0 = 0.05
inv.f = (num.persons - 6*2)*(1-alpha.0)/(6*(num.persons-6-1))
c(output,qf(inv.f, df1 = 6, df2 = num.persons - 6*2))
