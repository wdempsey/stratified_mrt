source("./functions.R")

N =  c(2,2)
num.persons = 15
T = 720
window.length = 60
P.treat = P.0 = P = matrix(c(0.6,0.4,0.04,0.96), byrow = TRUE, nrow = 2, ncol = 2)
eig.P = eigen(P)
pi = (eig.P$vectors%*%diag(c(1,0))%*%solve(eig.P$vectors))[1,]
tau = c(0.6,0.6)
lambda = 0.3
min.p = 0.001
max.p = 0.999
days = num.days = 10
percent.Delta = 0.1
beta = pi[1]*(1-percent.Delta)