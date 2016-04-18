### Estimation procedure given the dataset
people = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

colnames(people) = c("person", "day","t","Y.t","A.t","X.t","rho.t","I.t")

Y.t.person = people[,4]
A.t.person = people[,5]
X.t.person = people[,6]
rho.t.person = people[,7]

B.t.person = t(Vectorize(cov.gen)(people[,2]*people[,3]))

rho = 0.9

Z.t.person = B.t.person*matrix(rep(people[,4]-rho,3), byrow=TRUE, ncol = 3)

cov.t.person = cbind(B.t.person,Z.t.person)

log.weights = A.t.person*(log(rho) - log(rho.t.person)) + (1-A.t.person)*(log(1-rho) - log(1-rho.t.person))

fit.people = lm(Y.t.person~cov.t.person:as.factor(X.t.person)-1,weights = exp(log.weights))
fit.people2 = lm(Y.t.person~B.t.person:as.factor(X.t.person)-1,weights = exp(log.weights))

Covariates = model.matrix(fit.people)

XWX = foreach(i=1:num.persons, .combine = "+") %do% extract.tXWX(Covariates,people,log.weights,i)

Middle = foreach(person=1:num.persons, .combine = "+") %do% M.function(Covariates, people, log.weights, person,XWX,fit.people)

Sigma = solve(XWX,Middle)%*%solve(XWX)

(fit.people$coefficients[4:6]%*%solve(Sigma[4:6,4:6], fit.people$coefficients[4:6])+fit.people$coefficients[10:12]%*%solve(Sigma[10:12,10:12], fit.people$coefficients[10:12]))/num.persons

alpha.0 = 0.05
inv.f = (num.persons - 6*2)*(1-alpha.0)/(6*(num.persons-6-1))
qf(inv.f, df1 = 6, df2 = num.persons - 6*2)

# # W.hat = foreach(t=1:length(Y.t.person), .combine = "+") %do% fit.people$residuals[t]*exp(log.weights)[t]*Z.t.person[t,]
# # Q.hat = foreach(t=1:length(Y.t.person), .combine = "+") %do% outer(Z.t.person[t,],Z.t.person[t,])
# # 
# # 
# # Sigma = solve(Q.hat,outer(W.hat,W.hat))%*%solve(Q.hat)
# # 
# beta1 = fit.people$coefficients[c(5:7)]
# Sigma1 = vcov(fit.people)[c(4:6),c(4:6)]
# beta2 = fit.people$coefficients[c(12:14)]
# Sigma2 = vcov(fit.people)[c(10:12),c(10:12)]
# 
# 
# beta1%*%solve(Sigma1,beta1)+beta2%*%solve(Sigma2,beta2)
# num.persons*qf((num.persons-2*(3+3))*(1-0.05)/(2*3*(num.persons-2*3-1)),df1=2*3,df2=num.persons-2*(3+3))