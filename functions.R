### Functions associated with Sample Size Calculations
require(foreach); require(TTR)

### Simulation functions
rand.probs <- function(X.t, H.t, T, N, pi, tau, lambda, min.p, max.p) {
  ## Calculate randomization probabilities given pi, x.t, lambda, 
  ## T, and N 
  power = length(H.t$X):1
  rho.t = max(min((N[X.t] - sum(((1-lambda^power)*H.t$rho + lambda^power*H.t$A)*(H.t$X == X.t & H.t$I == 1)))/ (1 + (T-(max(power)+1))*pi[X.t]*tau[X.t]),max.p),min.p)
  return(rho.t)
}

calc.Ptreat <- function(pi, P, effect) {
  ## Calculate the Markov chain under treatment given pi, P, and effect
  pi2 = c(max(min(pi[1]-effect,1),0), max(min(pi[2]-effect,1),0))
  
  eig.P2 = eigen(rbind(pi2,pi2))
  
  test <- function(x) {
    sum((P - eig.P2$vectors%*%diag(c(1,x))%*%solve(eig.P2$vectors))^2)
  }
  
  P.treat = eig.P2$vectors%*%diag(c(1,optim(0.5, test,lower=0,upper = 1)$par))%*%solve(eig.P2$vectors)
  return(P.treat) 
}

daily.sim <- function(N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p) {
  ## Simulate the Markov chain given length (T), transition matrix (P), 
  ## and initial point (init)
  X.t = I.t = A.t = rho.t = vector(length = T)
  X.t[1] = sample(1:length(pi), size = 1, prob=pi)
  I.t[1] = rbinom(n=1,size = 1, prob=tau[X.t])
  if(I.t[1] == 1) {
    rho.t[1] = max(min(N[X.t]/((1 + (T - 1)*pi[X.t]*tau[X.t])),max.p),min.p)  
  } else{rho.t[1] = 0}
  A.t[1] = rbinom(n=1,size=1, prob = rho.t)
  H.t = list("X"=X.t[1],"A" = A.t[1], "I" = I.t[1], "rho" =rho.t[1])
  t = 2
  while (t <= T) {
    if(A.t[t-1] ==0) {
      X.t[t] = sample(1:nrow(P.0), size = 1, prob = P.0[X.t[t-1],])  
      I.t[t] = rbinom(n=1,size = 1, prob=tau[X.t])
      if( I.t[t] == 1) {
        rho.t[t] = rand.probs(X.t[t], H.t, T, N, pi, tau, lambda, min.p, max.p)
      } else ( rho.t[t] = 0 )
      A.t[t] = rbinom(n=1,size=1, prob = rho.t[t])
      H.t = list("X"=X.t[1:t],"A" = A.t[1:t], "I" = I.t[1:t], "rho" =rho.t[1:t])
      t = t+1
    } else {
      for(t.prime in t:(t+window.length-1)){
        X.t[t.prime] = sample(1:nrow(P.0), size = 1, prob = P.treat[X.t[t.prime-1],])  
      }
      I.t[t:(t+window.length-1)] = 0
      rho.t[t:(t+window.length-1)] = 0
      A.t[t:(t+window.length-1)] = 0
      H.t = list("X"=X.t[1:(t+window.length-1)],"A" = A.t[1:(t+window.length-1)], 
                 "I" = I.t[1:(t+window.length-1)], "rho" =rho.t[1:(t+window.length-1)])
      t = t+window.length
    }
  }
  return(H.t)
}

daily.sim_c <- compiler::cmpfun(daily.sim)



daily.data <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p){
  # Generate the daily data for a participant given all the inputs!
  inside.fn <- function(day) {
    P.treat = calc.Ptreat(pi,P,daily.treat[day])
    H.t = daily.sim(N, pi, tau, P.0, P.treat, T+window.length, window.length, min.p, max.p)
    Y.t = SMA(H.t$X==1,window.length); Y.t = Y.t[(window.length+1):length(Y.t)]
    data = cbind(day,1:T,Y.t,H.t$A[1:T],H.t$X[1:T], H.t$rho[1:T],H.t$I[1:T])
    return(data[data[,7] == 1,]  )
  }
  return(inside.fn)
}

full.trial.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  foreach(i=1:length(daily.treat), .combine = "rbind") %do% daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)(i)
}

MRT.sim <- function(num.people, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  # Do the trial across people!!
  foreach(i=1:num.people, .combine = "rbind") %do% cbind(i,full.trial.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p))
#   foreach(i=1:num.people, .combine = "rbind") %do% full.trial.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) # no cbind saves 0.32 seconds
}

cov.gen <-  function(t) {
  # For each t generate Z.t
  cov.t = c(1, floor((t-1)/T), floor((t-1)/T)^2)
  return(cov.t)
}

### Sample size calculations

find.d <- function(bar.d, init.d, max.d, Z.t, num.days) {
  # Find the quadratic terms given inputs
  D.star = num.days -1
  d = vector(length = 3)
  
  d[1] = 0
  d[3] = D.star * bar.d * solve(D.star^2 * ( D.star^2/3 - max.d))
  d[2] = -2 * d[3] * max.d
  
  ### Fix the scaling to get right average
  d = (ncol(Z.t) * bar.d / sum(d%*%Z.t)) * d
  
  return(d)
}

rho.function <- function(x, N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p) {
  ## Return the randomization probability for one simulated day
  H.t = daily.sim(N, pi, tau, P.0, P.treat, T+window.length, window.length, min.p, max.p)
  return(H.t$rho*(H.t$I == 1)*(H.t$X == x))
}

var.function <- function(x, N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p) {
  ## Return the variance of the randomization probability for one simulated day
  H.t = daily.sim(N, pi, tau, P.0, P.treat, T+window.length, window.length, min.p, max.p)
  Y.t = SMA(H.t$X==1,window.length); Y.t = Y.t[(window.length+1):length(Y.t)]
  temp = Y.t*(H.t$I[1:T] == 1)*(H.t$X[1:T] == x)
  return(var(temp[temp!=0]))
}

sample.size <- function(ss.param,p,q,alpha.0 = 0.05,beta.0 = 0.2, max.iters = 10000) {
  ## Sample size calculation given the ss.param, p, and q
  
  p.calc <- function(N) {
    c.N = N*ss.param
    
    # inv.q  = (N-q-p)*(1-alpha.0)/ (p*(N-q-1))
    df1 = p
    df2 = N-q-p
    
    inv.f = qf(1-alpha.0,df1,df2)
    
    return((p*(N-q-1))/(N-q-p)*pf(inv.f,df1,df2,ncp = c.N))
  }
  
  
  max.iters = 10000
  N = 100
  i = 1
  
  while (i < max.iters ) {
    if ( p.calc(N) < 1- beta.0) {
      if (p.calc(N-1) > 1 - beta.0) {
        break
      } else {N = N - 1}
    } else if (p.calc(N) > 1 - beta.0) {
      N = N+1
    } else if (p.calc(N) == 1 - beta.0) {break}
    i = i+1
  }
  
  return(N)
}

### Estimation

extract.tXWX <- function(cov, data, log.weights, person) {
  X = cov[data[,1]==person,]
  W = exp(log.weights[data[,1]==person])
  t(X)%*%diag(W)%*%X
}


M.function <- function(cov, data, log.weights, person, XWX, fit) {
  X = cov[data[,1]==person,]
  W = exp(log.weights[data[,1]==person])
  H = X%*%solve(XWX, t(X)%*%diag(W))
  e = fit$residuals[data[,1]==person]
  Q.test = solve(diag(nrow(H)) - H,e)
  M.i = t(X)%*%diag(W)%*%outer(Q.test,Q.test)%*%diag(W)%*%X
  return(M.i)
}

estimation <- function(people) {
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
  
  Covariates = model.matrix(fit.people)
  
  XWX = foreach(i=1:num.persons, .combine = "+") %do% extract.tXWX(Covariates,people,log.weights,i)
  
  Middle = foreach(person=1:num.persons, .combine = "+") %do% M.function(Covariates, people, log.weights, person,XWX,fit.people)
  
  Sigma = solve(XWX,Middle)%*%solve(XWX)
  
  entries = c(4:6,10:12)
  output = (fit.people$coefficients[entries]%*%solve(Sigma[entries,entries], fit.people$coefficients[entries]))/num.persons
  return(output)
}

estimation.simulation <- function(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  
  people = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)
  
  output = estimation(people)
  
  alpha.0 = 0.05
  inv.f = (num.persons - 6*2)*(1-alpha.0)/(6*(num.persons-6-1))
  return ( output > qf(inv.f, df1 = 6, df2 = num.persons - 6*2))
}
