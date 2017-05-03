### Functions associated with Sample Size Calculations
require(foreach); require(TTR); require(zoo); require(expm)

### Simulation functions
rand.probs <- function(X.t, H.t, T, N, pi, tau, lambda, min.p, max.p) {
  ## Calculate randomization probabilities given pi, x.t, lambda,
  ## T, and N
  power = length(H.t$X):1
  remaining.time = T-(max(power)+1)
  if(remaining.time < 0) {
    true.rem.time = 0
  } else if(remaining.time - N[X.t]*60 < 60) {
    true.rem.time = remaining.time*pi[X.t]*tau[X.t]
  } else if(remaining.time - N[X.t]*60 < 120) {
    true.rem.time = (remaining.time-60)*pi[X.t]*tau[X.t]
  } else {
    true.rem.time = (remaining.time - 120)*pi[X.t]*tau[X.t]
  }
  #   rho.t = max(min((N[X.t] - sum(((1-lambda^power)*H.t$rho + lambda^power*H.t$A)*(H.t$X == X.t & H.t$I == 1)))/ (1 + (T-(max(power)+1))*pi[X.t]*tau[X.t]),max.p),min.p)
  rho.t = max(min((N[X.t] - sum(((1-lambda^power)*H.t$rho + lambda^power*H.t$A)*(H.t$X == X.t & H.t$I == 1)))/ (1 + true.rem.time),max.p),min.p)
  return(rho.t)
}

potential.effects <- function(P,window.length) {
  max.delta = 1-diag(P)
  min.delta = -diag(P)
  
  delta.1.range = seq(min.delta[1],max.delta[1],0.005)
  delta.2.range = seq(min.delta[2],max.delta[2],0.005)
  
  direct.treat.X1 = direct.treat.X2 = matrix(nrow = length(delta.1.range), ncol = length(delta.2.range))
  
  for (i in 1:length(delta.1.range)) {
    for(j in 1:length(delta.2.range)) {
      
      delta.1 = delta.1.range[i]; delta.2 = delta.2.range[j]
      
      P.treat = P + matrix(c(delta.1,-delta.1,-delta.2,delta.2), nrow = 2, ncol = 2, byrow = TRUE)
      
      direct.result = list()
      direct.result$A1 = direct.result$A0 = 0
      
      for (k in 1:window.length) {
        direct.result$A0 = direct.result$A0 + (P%^%k)[,2]
        direct.result$A1 = direct.result$A1 + (P.treat%^%k)[,2]
      }
      
      temp = (direct.result$A1 -
                direct.result$A0)/window.length
      
      direct.treat.X1[i,j] = temp[1]; direct.treat.X2[i,j] = temp[2]
      
    }
  }
  
  data.frame.treat = matrix(nrow = length(delta.1.range)*length(delta.2.range), ncol = 4)
  
  for (i in 1:length(delta.1.range)) {
    for(j in 1:length(delta.2.range)) {
      data.frame.treat[(i-1)*length(delta.2.range)+j,] = c(delta.1.range[i],delta.2.range[j],direct.treat.X1[i,j], direct.treat.X2[i,j])
    }
  }
  
  data.frame.treat = data.frame(data.frame.treat)
  
  names(data.frame.treat) = c("delta.1", "delta.2", "treat.X1", "treat.X2")
  
  return(data.frame.treat)
}

calc.Ptreat <- function(P, effect, treatment.data, tol) {
  
  if (all(effect == 0) == TRUE) {
    return(P)
  } else {
    
    obs1 = which(treatment.data$treat.X1 > effect[1] - tol & treatment.data$treat.X1 < effect[1] + tol)
    obs2 = which(treatment.data$treat.X2 > effect[2] - tol & treatment.data$treat.X2 < effect[2] + tol)
    
    if(length(obs1) == 0 | length(obs2) == 0) {stop("no observations of this effect combination with set tolerance level")}
    
    fit.obs1 = lm(treatment.data$delta.2[obs1]~treatment.data$delta.1[obs1])
    fit.obs2 = lm(treatment.data$delta.2[obs2]~treatment.data$delta.1[obs2])
    
    delta.1.intersect = -(fit.obs1$coefficients[1] - fit.obs2$coefficients[1])/(fit.obs1$coefficients[2] - fit.obs2$coefficients[2])
    delta.2.intersect = fit.obs1$coefficients[1] + fit.obs1$coefficients[2]*delta.1.intersect
    
    P.treat = P + matrix(c(delta.1.intersect,-delta.1.intersect,-delta.2.intersect,delta.2.intersect), nrow = 2, ncol = 2, byrow = TRUE)
    
    ## Quick check on treatment effect
    
    direct.result = list()
    direct.result$A1 = direct.result$A0 = 0
    
    for (k in 1:60) {
      direct.result$A0 = direct.result$A0 + (P%^%k)[,2]
      direct.result$A1 = direct.result$A1 + (P.treat%^%k)[,2]
    }
    
    temp = (direct.result$A1 - direct.result$A0)/60
    
    if(!all(abs(temp - effect)< tol)){
      stop("no observations of this effect combination with set tolerance level")
    }
    
    
    return(P.treat)
  }
}

daily.sim <- function(N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p) {
  ## Simulate the Markov chain given length (T), transition matrix (P),
  ## and initial point (init)
  X.t = I.t = A.t = rho.t = vector(length = T)
  X.t[1] = sample(1:length(pi), size = 1, prob=pi)
  I.t[1] = rbinom(n=1,size = 1, prob=tau[X.t])
  if(I.t[1] == 1) {
    rho.t[1] = max(min(N[X.t[1]]/((1 + (T - 120 - 1)*pi[X.t]*tau[X.t[1]])),max.p),min.p)
  } else{rho.t[1] = 0}
  A.t[1] = rbinom(n=1,size=1, prob = rho.t)
  H.t = list("X"=X.t[1],"A" = A.t[1], "I" = I.t[1], "rho" =rho.t[1])
  t = 2
  while (t <= T+window.length) {
    if(A.t[t-1] ==0) {
      X.t[t] = sample(1:nrow(P.0), size = 1, prob = P.0[X.t[t-1],])
      if(t > T) {
        I.t[t] = 0; rho.t[t] = 0
      } else{
        I.t[t] = rbinom(n=1,size = 1, prob=tau[X.t[t]])
        if( I.t[t] == 1) {
          rho.t[t] = rand.probs(X.t[t], H.t, T, N, pi, tau, lambda, min.p, max.p)
        } else ( rho.t[t] = 0 )
      }
      A.t[t] = rbinom(n=1,size=1, prob = rho.t[t])
      H.t = list("X"=X.t[1:t],"A" = A.t[1:t], "I" = I.t[1:t], "rho" =rho.t[1:t])
      t = t+1
    } else {
      I.t[t:(t+window.length-1)] = 0
      rho.t[t:(t+window.length-1)] = 0
      A.t[t:(t+window.length-1)] = 0
      for(t.prime in t:(t+window.length-1)){
        X.t[t.prime] = sample(1:nrow(P.0), size = 1, prob = P.treat[X.t[t.prime-1],])
      }
      H.t = list("X"=X.t[1:(t+window.length-1)],"A" = A.t[1:(t+window.length-1)],
                 "I" = I.t[1:(t+window.length-1)], "rho" =rho.t[1:(t+window.length-1)])
      t = t+window.length
    }
  }
  return(H.t)
}

daily.data <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data){
  # Generate the daily data for a participant given all the inputs!
  
  inside.fn <- function(day) {
    effect = rep(daily.treat[day],2)
    P.treat = calc.Ptreat(P.0,effect,treatment.data, tol=10^(-2))
    H.t = daily.sim(N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
    Y.t = SMA(H.t$X==2,window.length); Y.t = Y.t[(window.length+1):(length(Y.t))]
    #Y.t = sim1_Y(H.t,day,d)[1:T]
    prob.gamma = rollapply(1-H.t$rho, window.length, FUN = prod); prob.gamma = prob.gamma[-1]
    prob.nu = rollapply((H.t$A==0),window.length, FUN = prod); prob.nu = prob.nu[-1]
    psi.t = prob.nu/prob.gamma
    data = cbind(day,1:T,Y.t,H.t$A[1:T],H.t$X[1:T], H.t$rho[1:T],H.t$I[1:T], psi.t)
    return(data[data[,7] == 1 & data[,8] > 0,])
  }
  return(inside.fn)
}

full.trial.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  foreach(i=1:length(daily.treat), .combine = "rbind", .packages = c("foreach", "TTR","expm","zoo")) %dopar% daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)(i)
}

MRT.sim <- function(num.people, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data) {
  # Do the trial across people!!
  output = foreach(i=1:num.people, .combine = "rbind", .packages = c("foreach", "TTR","expm","zoo")) %dopar% cbind(i,full.trial.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data))
  colnames(output) = c("person", "day", "t", "Y.t","A.t","X,t", "rho.t", "I.t","psi.t")
  return(output)
}

f.t <-  function(t,X.t) {
  # For each t generate X.t
  cov.t = c(1, floor((t-1)/T), floor((t-1)/T)^2)
  return(rep(cov.t,2)*c(rep(X.t==1,3),rep(X.t==2,3)))
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
  
  Y.t.person = people[,4]
  A.t.person = people[,5]
  X.t.person = people[,6]
  rho.t.person = people[,7]
  psi.t.person = people[,9]
  
  B.t.person = t(Vectorize(cov.gen)((people[,2]-1)*T + people[,3]))
  
  # Set of possible weights depending on unique X.t
  set.rho = foreach(lvl=1:length(unique(X.t.person)), .combine = "c", .packages = c("foreach", "TTR","expm","zoo")) %dopar% mean(rho.t.person[X.t.person==lvl], na.rm = TRUE)
  
  rho = unlist(lapply(X.t.person,tilde.p))
  
  Z.t.person = B.t.person*matrix(rep(A.t.person-rho,3), ncol = 3)
  
  cov.t.person = cbind(B.t.person,Z.t.person)
  
  log.weights = A.t.person*(log(rho) - log(rho.t.person)) + (1-A.t.person)*(log(1-rho) - log(1-rho.t.person)) + log(psi.t.person)
  
  fit.people = lm(Y.t.person~(B.t.person+Z.t.person):as.factor(X.t.person)-1,weights = exp(log.weights))
  
  Covariates = model.matrix(fit.people)
  
  num.persons = length(unique(people[,1]))
  
  XWX = foreach(i=1:num.persons, .combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dopar% extract.tXWX(Covariates,people,log.weights,i)
  
  Middle = foreach(person=1:num.persons, .combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dopar% M.function(Covariates, people, log.weights, person,XWX,fit.people)
  
  entries1 = c(7:9)
  entries2 = c(10:12)
  entries = c(entries1,entries2)
  
  Sigma = solve(XWX,Middle)%*%solve(XWX)
  
  output = (fit.people$coefficients[entries]%*%solve(Sigma[entries,entries], fit.people$coefficients[entries]))
  
  return(output)
}

estimation.simulation <- function(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data) {
  
  people = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data)
  
  output = estimation(people)
  
  alpha.0 = 0.05; p = 6; q = 6;
  
  multiple = p*(num.persons-q-1)/(num.persons-p-q)
  
  return (output>multiple*qf((1-alpha.0), df1 = p, df2 = num.persons - p - q))
}

#### SS-Calculation functions
tilde.p <- function(X.t) {
  ## Constant fn of t, stratification across X.t 
  # N[X.t]/((T-60*N[X.t])*pi[X.t]*tau[X.t])
  ## Constant fn of t, and X.t
  N[1]/((T-60*N[1])*pi[1]*tau[1])*pi[1]+N[2]/((T-60*N[2])*pi[2]*tau[2])*pi[2]
}