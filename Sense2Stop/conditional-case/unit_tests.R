# Check for each time treated,  you are unavailable for next 60 minutes
set.seed(19731)
temp = daily.sim(N, pi, P, P.treat.list[[5]], T, window.length, min.p, max.p) 
treated_obs = which(temp$A == 1)
for (obs in treated_obs) {
  print(any(temp$I[(obs+1):(obs+window.length)] == 1))
}

# Set P to be 2->3, 3->3, and P.treat to 3->2, 2->2,
# To make sure the transitions are when they should be
set.seed(19731)
P.test = P; 
P.test[3,3] = 0; P.test[3,] = P.test[3,]/sum(P.test[3,])
P.test[6,6] = 0; P.test[6,] = P.test[6,]/sum(P.test[6,])
eig.P.test = eigen(P.test)
pi.test = (eig.P$vectors%*%diag(c(1,rep(0,5)))%*%solve(eig.P$vectors))[1,]  # Stationary distribution for P
P.treat.test = P.treat;
P.treat.test[3,3] = 1; P.treat.test[3,1] = 0; P.treat.test[3,4] = 0
P.treat.test[6,6] = 1; P.treat.test[6,1] = 0; P.treat.test[6,4] = 0
temp = daily.sim(N, pi.test, P.test, P.treat.test, T, window.length, min.p, max.p) 
treated_obs = which(temp$A == 1)
# We should see state 5 or 2, then 60 times in state 6 or 3 respectively, then immediately change to 
# state 4 or 1 (didn't want to make it irreducible)
for (obs in treated_obs) {
  print(temp$X[(obs):(obs+window.length+1)])
}



num.people=100
temp = daily.data(N, pi, P, P.treat.list, T, window.length, min.p, max.p)(5)
for(iter in 1:1000) {
temp = rbind(temp, daily.data(N, pi, P.0, P.treat.list, T, window.length, min.p, max.p)(5))
}
output = as.data.frame(temp)
names(output) = c("day", "t", "Y.t", "A.t", "X.t", "I.t", "psi.t")
mean(output$Y.t[output$X.t == 2 & output$A.t == 1])
mean(output$Y.t[output$X.t == 2 & output$A.t == 0])


mean(output$Y.t[output$X.t == 5 & output$A.t == 1])
mean(output$Y.t[output$X.t == 5 & output$A.t == 0])

