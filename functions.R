### Functions associated with Sample Size Calculations

### Sample sequence of states

daily.sim <- function(init, P, T) {
  ## Simulate the Markov chain given length (T), transition matrix (P), 
  ## and initial point (init)
  X.t = vector(length = T)
  X.t[1] = init
  for(i in 2:T) {
    X.t[i] = sample(1:nrow(P), size = 1, prob = P[X.t[i-1],])
  }
  return(X.t)
  
}

rand.probs <- function(X.t, T, N, pi, lambda, min.p, max.p) {
  ## Calculate randomization probabilities given pi, x.t, lambda, 
  ## T, and N
  
  rho.t = matrix(nrow = length(pi), ncol = T)  
  A.t = matrix(0,nrow = length(pi), ncol = T)  
  rho.t[,1] = N/(T*pi)
  A.t[,1] = rbinom(length(pi), size = 1, prob = rho.t[,1])
  for(i in 2:T) {
    for(k in 1:length(pi)){
      if (N[k] == 0) {
        rho.t[k,i] = 0; A.t[k,i] = 0
      } else {
        rho.t[k,i] = max(min((N[k] - sum((lambda*rho.t[k,1:T<i]+(1-lambda)*A.t[k,1:T<i])*(X.t[1:T<i] == k)))/ ((T-i+1)*pi[k]),max.p),min.p)
        A.t[k,i] = rbinom(1, size = 1, prob = rho.t[k,i])
      }
    }
  }
  return(rho.t)
}

inverse.expectation <- function(num.iters, pi, P,T, N, min.p, max.p) {

  inv.var.est = matrix(0,nrow = 2, ncol = T)
  counts = matrix(0,nrow = 2, ncol = T)
  
  for(i in 1:num.iters) {
    init = sample(1:3,size = 1, prob = pi)
    X.t_i = daily.sim(init, P, T)
    rho.t_i = rand.probs(X.t_i, T, N, pi, 1, min.p, max.p)
    for(k in 1:2) {
      inv.var.est[k,X.t_i == k] = inv.var.est[k,X.t_i == k]  + 1/(rho.t_i[k,X.t_i == k]*(1-rho.t_i[k,X.t_i == k]))
      counts[k,X.t_i ==k] = counts[k,X.t_i ==k] + 1
    }
  }
  
  return(inv.var.est/counts)
}

  
### Conditional sample size calculation
# We are averaging over E[ (rho_t*(1-rho_t))^{-1} | X_t = x] for now

sample.size <- function(gamma, inv.var, avgtreatment,T,days) {
  
  # Number of levels for time-varying covariate
  k = 2
  
  # Observation times
  t = seq(1,days*T)
  
  # For each t generate Z.t
  cov.gen <-  function(t) {
    cov.t = c(1, floor((t-1)/T), floor((t-1)/T)^2)
    return(cov.t)
  }
  
  Z.t = Vectorize(cov.gen)(t)
  
  # Standardized treatment effect (constant over levels)
  bar.d = rep(avgtreatment,2)
  init.d = rep(0,2)
  max.d = 14/2
  
  c.star = max.d*2
  d = vector(length = 3)
  
  d[2] = max(t)*bar.d[1]/sum(c(-1/c.star, 1, (1+1/c.star))%*%Z.t)
  d[1] = d[2]*(1+1/c.star)
  d[3] = d[2]*(-1/c.star)
  
  ### Fix the scaling to get right average
  d = (length(t) * bar.d[1] / sum(d%*%Z.t)) * d
  
  # Non-central chi-squared parameters
  
  W = Q = list(); inv.list = list()
  
  for (i in 1:k) {
    inv.list[[i]] = rep(inv.var.est[i,],days)
    W[[i]] = Q[[i]]= matrix(0,nrow = 3, ncol = 3) 
  }
  
  for (i in 1:length(t)) {
    for (j in 1:k) {
      Q[[j]] = Q[[j]] + (Z.t[,i] %*% t(Z.t[,i]))
      W[[j]] = W[[j]] + (Z.t[,i] %*% t(Z.t[,i]))*(inv.list[[j]][i])
    }
  }
  
  
  p.calc <- function(N) {
    c.N = 0
    
    for (j in 1:k) {
      tot = Q[[j]]%*%solve(W[[j]],Q[[j]])
      c.N = c.N + N*(t(d)%*%tot%*%d)*gamma^2
    }
    
    # Calc the F-distribution
    
    q = 3
    p = 3
    alpha.0 = 0.05
    
    # inv.q  = (N-q-p)*(1-alpha.0)/ (p*(N-q-1))
    df1 = p
    df2 = N-q-p
    
    inv.f = qf(1-alpha.0,df1,df2)
    
    return((p*(N-q-1))/(N-q-p)*pf(inv.f,df1,df2,ncp = c.N))
  }
  
  
  beta.0 = 0.2
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

find.d <- function(bar.d, init.d, max.d) {
  
  c.star = max.d*2
  d = vector(length = 3)
  
  d[2] = max(t)*bar.d[1]/sum(c(-1/c.star, 1, (1+1/c.star))%*%Z.t)
  d[1] = d[2]*(1+1/c.star)
  d[3] = d[2]*(-1/c.star)
  
  ### Fix the scaling to get right average
  d = (length(t) * bar.d[1] / sum(d%*%Z.t)) * d
  
  return(d)
}