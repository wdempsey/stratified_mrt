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
std.beta = 0.0025/0.004

# Observation times
t = seq(1,days*T)

Z.t = Vectorize(cov.gen)(t)

d = find.d(bar.d = std.beta*0.004,init.d = 0,max.d = 5,Z.t = Z.t,num.days=num.days)

t.daily = seq(1,days*T,T)

Z.t.daily = Vectorize(cov.gen)(t.daily)

daily.treat = t(Z.t.daily)%*%d

par(mar = c(4,4,1,1)+0.1, mfrow = c(1,1))
plot(1:days,daily.treat, xlab = "Day", ylab = "Non-standardized daily effect",axes = FALSE, xlim = c(0,10))
abline(v = 6, lty = 2)
axis(side = 1); axis(side = 2)

daily.treat[6]

calc.Ptreat(pi,P,daily.treat[6])
