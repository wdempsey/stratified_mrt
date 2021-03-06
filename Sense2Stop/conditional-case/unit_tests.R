source('./setup.R'); source("./functions.R")
# Make a no treatment effect setting
P.treat.list = list()
for(day in 1:num.days) {
  P.treat.list[[day]] = P
}

# Check 1: Each time treated,  you are unavailable for next 60 minutes
set.seed(19731)
temp = daily.sim(N, pi, P, P.treat.list[[5]], T, window.length, min.p, max.p) 
treated_obs = which(temp$A == 1)
# Should report all false for availability
for (obs in treated_obs) {
  print(any(temp$I[(obs+1):(obs+window.length)] == 1))
}

# Check 2: Check if transitions are generated properly
# Set P to be 2->3, 3->3, and P.treat to 3->2, 2->2,
# To make sure the transitions are when they should be
set.seed(19731)
P.test = P; 
P.test[3,3] = 0; P.test[3,] = P.test[3,]/sum(P.test[3,])
P.test[6,6] = 0; P.test[6,] = P.test[6,]/sum(P.test[6,])
eig.P.test = eigen(P.test)
pi.test = (eig.P$vectors%*%diag(c(1,rep(0,5)))%*%solve(eig.P$vectors))[1,]  # Stationary distribution for P
P.treat.test = P.treat.list[[5]];
P.treat.test[3,3] = 1; P.treat.test[3,1] = 0; P.treat.test[3,4] = 0
P.treat.test[6,6] = 1; P.treat.test[6,1] = 0; P.treat.test[6,4] = 0
temp = daily.sim(N, pi.test, P.test, P.treat.test, T, window.length, min.p, max.p) 
treated_obs = which(temp$A == 1)
Y.t = SMA(is.element(temp$X,c(4,5,6)),window.length); Y.t = Y.t[(window.length+1):(length(Y.t))]
# We should see state 5 or 2, then 60 times in state 6 or 3 respectively, then immediately change to 
# state 4 or 1 (didn't want to make it irreducible)
for (obs in treated_obs) {
  print(temp$X[(obs):(obs+window.length+1)])
  print(Y.t[obs])
}


# Check 3: Check Y.t computed properly
set.seed(19731)
temp = daily.sim(N, pi, P, P.treat.list[[5]], T, window.length, min.p, max.p) 
Y.t.forloop = vector(length = T)
for(t in 1:T) {
  Y.t.forloop[t] = mean(is.element(temp$X[(t+1):(t+window.length)],c(4:6)))
}
Y.t.sma = SMA(is.element(temp$X,c(4,5,6)),window.length); Y.t.sma = Y.t.sma[(window.length+1):(length(Y.t.sma))]
all((Y.t.forloop-Y.t.sma)==0) # Shows that SMA works properly

# Check 4: Rollapply is computing pr(no action in next 60 minutes) properly
# and indicator of no action in next 60 minutes
prob.gamma.forloop = prob.nu.forloop = vector(length = T)
for(t in 1:T) {
  prob.gamma.forloop[t] = prod(1-temp$rho[(t+1):(t+window.length-1)])
  prob.nu.forloop[t] = prod(temp$A[(t+1):(t+window.length-1)]==0)
}
prob.gamma.rollapply = rollapply(1-temp$rho, window.length-1, FUN = prod); prob.gamma.rollapply = prob.gamma.rollapply[-c(1, length(prob.gamma.rollapply))]
all((prob.gamma.forloop - prob.gamma.rollapply)==0) # Shows that rollapply works properly
prob.nu.rollapply = rollapply((temp$A==0),window.length-1, FUN = prod); prob.nu.rollapply = prob.nu.rollapply[-c(1, length(prob.nu.rollapply))]
all((prob.nu.forloop - prob.nu.rollapply)==0) # Shows that rollapply works properly


## Check 5: Check Y.t under no effect is near truth ##
# Compute true means
total = 0.0
for (k in 1:window.length) {
  total = total + (rowSums((P%^%k)[c(2,5),4:6]))
}

set.seed(19371)
total.iter = 1000
# Run 5000 person-days under null hypothesis
res1 = res2 = res3 = res4 = rep(0,0)
for (iter in 1:total.iter) {
  mod.iter = iter%%10
  if(mod.iter == 0) {mod.iter = 10}
  temp.person = daily.sim(N, pi, P, P.treat.list[[mod.iter]], T, window.length, min.p, max.p)
  Y.t = SMA(is.element(temp.person$X,c(4,5,6)),window.length); 
  Y.t = Y.t[(window.length+1):(length(Y.t))]
  res1 = c(res1, Y.t[temp.person$A == 1 & temp.person$X == 2 & temp.person$I == 1])
  res2 = c(res2, Y.t[temp.person$A == 0 & temp.person$X == 2 & temp.person$I == 1])
  
  res3 = c(res3, Y.t[temp.person$A == 1 & temp.person$X == 5 & temp.person$I == 1])
  res4 = c(res4, Y.t[temp.person$A == 0 & temp.person$X == 5 & temp.person$I == 1])
}

# Out of 1000 days, how many (roughly) in each block
length(res1); length(res2); length(res3); length(res4)
# As expected, we have roughly 1500/1000 = 1.5 treatments
# on average across days
print(c(mean(res1), mean(res2), total[1]/window.length))
print(c(mean(res3), mean(res4), total[2]/window.length))

## Check full.trial per person is running correctly
set.seed(19371)
total.people = 100
# Run 5000 person-days under null hypothesis
output = rep(0,0)
for (person in 1:total.people) {
  for (day in 1:10) {
    temp.person = cbind(person, daily.data(N, pi, P, P.treat.list, T, window.length, min.p, max.p)(day))
    output = rbind(output, temp.person)
  }
}
output = data.frame(output)
colnames(output) = c("person", "day", "t", "Y.t","A.t","X.t", "rho.t", "I.t","psi.t")
sum(output$A.t == 1 & output$X.t == 2 & output$I.t == 1 & output$psi.t > 0) # Should be approx 1.5k
sum(output$A.t == 0 & output$X.t == 2 & output$I.t == 1 & output$psi.t > 0)
sum(output$A.t == 1 & output$X.t == 5 & output$I.t == 1 & output$psi.t > 0) # Should be approx 1.5k
sum(output$A.t == 0 & output$X.t == 5 & output$I.t == 1 & output$psi.t > 0)

# I should see that the means are find when I don't control for psi.t
print(paste("True mean is", round(total[1]/window.length,4)))
mean(output$Y.t[output$A.t == 1 & output$X.t == 2])
mean(output$Y.t[output$A.t == 0 & output$X.t == 2])

print(paste("True mean is", round(total[2]/window.length,4)))
mean(output$Y.t[output$A.t == 1 & output$X.t == 5])
mean(output$Y.t[output$A.t == 0 & output$X.t == 5])

prod((1-output$rho.t[2:(2+window.length)]))
## Now, the data we analyze is only 
## when participant is available and no subsequent treatment
## that is I_t = 1, psi_t > 0.  But then we need to adjust the 
## expectations to see if unbiased.
## This will add variance to our estimators, but will be unbiased.
print(paste("True mean is", round(total[1]/window.length,4)))
mean((output$psi.t*output$Y.t)[output$A.t == 1 & output$X.t == 2])
mean((output$psi.t*output$Y.t)[output$A.t == 0 & output$X.t == 2])

print(paste("True mean is", round(total[2]/window.length,4)))
mean((output$psi.t*output$Y.t)[output$A.t == 1 & output$X.t == 5])
mean((output$psi.t*output$Y.t)[output$A.t == 0 & output$X.t == 5])

## SAME CHECK FOR MRT.SIM
## Check full.trial per person is running correctly
set.seed(19371)
total.people = 100
output = MRT.sim(num.people = total.people, N, pi, P, P.treat.list, T, window.length, min.p, max.p)
sum(output$A.t == 1 & output$X.t == 2 & output$I.t == 1 & output$psi.t > 0) # Should be approx 1.5k
sum(output$A.t == 0 & output$X.t == 2 & output$I.t == 1 & output$psi.t > 0)
sum(output$A.t == 1 & output$X.t == 5 & output$I.t == 1 & output$psi.t > 0) # Should be approx 1.5k
sum(output$A.t == 0 & output$X.t == 5 & output$I.t == 1 & output$psi.t > 0)

# I should see that the means are find when I don't control for psi.t
print(paste("True mean is", round(total[1]/window.length,4)))
mean(output$Y.t[output$A.t == 1 & output$X.t == 2])
mean(output$Y.t[output$A.t == 0 & output$X.t == 2])

print(paste("True mean is", round(total[2]/window.length,4)))
mean(output$Y.t[output$A.t == 1 & output$X.t == 5])
mean(output$Y.t[output$A.t == 0 & output$X.t == 5])

## Now, the data we analyze is only 
## when participant is available and no subsequent treatment
## that is I_t = 1, psi_t > 0.  But then we need to adjust the 
## expectations to see if unbiased.
## This will add variance to our estimators, but will be unbiased.
print(paste("True mean is", round(total[1]/window.length,4)))
sum((output$psi.t*output$Y.t)[output$A.t == 1 & output$X.t == 2])/sum(output$psi.t[output$A.t == 1 & output$X.t == 2])
sum((output$psi.t*output$Y.t)[output$A.t == 0 & output$X.t == 2])/sum(output$psi.t[output$A.t == 0 & output$X.t == 2])

print(paste("True mean is", round(total[2]/window.length,4)))
sum((output$psi.t*output$Y.t)[output$A.t == 1 & output$X.t == 5])/sum(output$psi.t[output$A.t == 1 & output$X.t == 5])
sum((output$psi.t*output$Y.t)[output$A.t == 0 & output$X.t == 5])/sum(output$psi.t[output$A.t == 0 & output$X.t == 5])

## Generate test statistic info under the null hypothesis
B.t.person = t(Vectorize(cov.gen)((output$day-1)*T + output$t))
