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

potential.effects <- function(P) {
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

      for (k in 1:60) {
        direct.result$A0 = direct.result$A0 + (P%^%k)[,1]
        direct.result$A1 = direct.result$A1 + (P.treat%^%k)[,1]
      }

      temp = (direct.result$A1 - direct.result$A0)/60

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
        direct.result$A0 = direct.result$A0 + (P%^%k)[,1]
        direct.result$A1 = direct.result$A1 + (P.treat%^%k)[,1]
    }

    temp = (direct.result$A1 - direct.result$A0)/60

    if(!all(abs(temp - effect)< tol)){
        stop("no observations of this effect combination with set tolerance level")
    }


  return(P.treat)
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

daily.sim_c <- compiler::cmpfun(daily.sim)


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

MRT.sim <- function(num.people, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  # Do the trial across people!!
  treatment.data = potential.effects(P)
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

# ss.parameters <- function(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
#     ## Return the randomization probability for one simulated day
#     full.sim = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)
#     full.sim = data.frame(full.sim); colnames(full.sim) = c("person","day", "t", "Y.t", "A.t", "X.t", "rho.t", "I.t", "psi.t")
#     full.sim$log.integrand = 2 * (-full.sim$A.t * log(full.sim$rho.t) -
#                                     (1-full.sim$A.t)*log(1-full.sim$rho.t)) +
#                                     log(full.sim$psi.t)
#     full.sim$integrand = exp(full.sim$log.integrand)
#     return(list("exp.est"=aggregate(cbind(integrand)~
#                                  X.t+day,
#                              data = full.sim,
#                              mean, na.rm = TRUE),
#                 "sigmasq"=var(full.sim$Y.t)))
# }

sample.size <- function(ss.param,p,q,alpha.0 = 0.05,beta.0 = 0.8, max.iters = 10000) {
  ## Sample size calculation given the ss.param, p, and q

  p.calc <- function(N) {
    c.N = N*ss.param

    # inv.q  = (N-q-p)*(1-alpha.0)/ (p*(N-q-1))
    df1 = p
    df2 = N-q-p

    inv.f = qf(1-alpha.0,df1,df2)

    return(pf(inv.f,df1,df2,ncp = c.N))
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

  Y.t.person = people[,4]
  A.t.person = people[,5]
  X.t.person = people[,6]
  rho.t.person = people[,7]
  psi.t.person = people[,9]

  B.t.person = t(Vectorize(cov.gen)((people[,2]-1)*T + people[,3]))

  # Set of possible weights depending on unique X.t
  set.rho = foreach(lvl=1:length(unique(X.t.person)), .combine = "c", .packages = c("foreach", "TTR","expm","zoo")) %dopar% mean(rho.t.person[X.t.person==lvl], na.rm = TRUE)

#   rho.val <- function(lvl) {
#     return(set.rho[lvl])
#   }

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

estimation.simulation <- function(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {

    people = MRT.sim(num.persons, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

    output = estimation(people)

    alpha.0 = 0.05; p = 6; q = 6;

    multiple = p*(num.persons-q-1)/(num.persons-p-1)

  return (output>multiple*qf((1-alpha.0), df1 = p, df2 = num.persons - p-q))
}


#### SS-Calculation functions
tilde.p <- function(X.t) {
  N[X.t]/((T-60*N[X.t])*pi[X.t]*tau[X.t])
}

gamma.st <- function(s,t) {
  ### Simple correlation model
  max(1-abs(s-t)*1.70*10^(-2),0)
}

ss.daily.data <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data){
  # Generate the daily data for a participant given all the inputs!

  inside.fn <- function(day) {
    effect = rep(daily.treat[day],2)
    P.treat = calc.Ptreat(P.0,effect,treatment.data, tol=10^(-2))
    bar.sigma = bar.sigmasq.fn(P,pi,P.treat,window.length)
    H.t = daily.sim(N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
    prob.gamma = rollapply(1-H.t$rho, window.length, FUN = prod); prob.gamma = prob.gamma[-1]
    prob.nu = rollapply((H.t$A==0),window.length, FUN = prod); prob.nu = prob.nu[-1]
    psi.t = prob.nu/prob.gamma
    data = cbind(day,1:(T+window.length),H.t$A,H.t$X, H.t$rho,H.t$I, c(psi.t,rep(0,window.length)))
    Q = W.first = W.second = 0
    obs.times = data[data[,6]==1 & data[,7] > 0,2]
    for(i in 1:(length(obs.times))) {
      ob.t = obs.times[i] ; f_at_obt = f.t(ob.t+ T*day, H.t$X[ob.t])
      p.t = tilde.p(H.t$X[ob.t])
      add.weight.t = (p.t/H.t$rho[ob.t])^(H.t$A[ob.t]) * ((1-p.t)/(1-H.t$rho[ob.t]))^(1-H.t$A[ob.t])
      sigmasq.ob.t = H.t$A[ob.t]*sigmasq.fn(P.treat,H.t$X[ob.t],window.length)+(1-H.t$A[ob.t])*sigmasq.fn(P,H.t$X[ob.t],window.length)
      Q = Q +
          psi.t[ob.t] * add.weight.t * (H.t$A[ob.t]-p.t)^2*outer(f_at_obt,f_at_obt)
      W.first = W.first +
          psi.t[ob.t]^2 * add.weight.t^2 * (H.t$A[ob.t]-p.t)^2*sigmasq.ob.t*outer(f_at_obt,f_at_obt)          
      diff = obs.times[-i]-ob.t
      other.times = obs.times[-i]
      close.times = other.times[abs(diff)<window.length]
      if(length(close.times) > 0) {
          for (j in 1:length(close.times)) {
              ob.s = close.times[j]; f_at_obs = f.t(ob.s+T*day,H.t$X[ob.s])
              p.s = tilde.p(H.t$X[ob.s])
              add.weight.s = (p.s/H.t$rho[ob.s])^(H.t$A[ob.s]) * ((1-p.s)/(1-H.t$rho[ob.s]))^(1-H.t$A[ob.s])
              W.second = W.second + bar.sigma *
                  (1-abs(ob.s-ob.t)/window.length)*
                  psi.t[ob.t] * add.weight.t * (H.t$A[ob.t]-p.t) *
                  psi.t[ob.s] * add.weight.s * (H.t$A[ob.s]-p.s) *
                  outer(f_at_obt,f_at_obs)
          }
      }
    }
    return(list("Q.day" = Q, "W.day" = W.first+W.second))
  }
  return(inside.fn)
}

full.trial.ss.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  W = Q = outer(f.t(1,1),f.t(1,1))*0
  for(i in 1:length(daily.treat)) {
    day.res = ss.daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)(i)
    W = W + day.res$W.day
    Q = Q + day.res$Q.day
  }
  return(rbind(Q,W))
}

ss.parameters <- function(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  treatment.data = potential.effects(P)
  output = foreach(i=1:num.iters,.combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dopar% full.trial.ss.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)/num.iters
  return(output)
}

sigmasq.fn <- function(P,current.state,window.length) {
    var.total = 0
    for (s in 1:window.length) {
        for (s.prime in 1:window.length) {
            delta = max(s,s.prime) - min(s,s.prime)
            var.total = var.total + (P%^%min(s,s.prime))[current.state,2] *
                (P%^%delta)[2,2] -
                (P%^%s)[current.state,2]*(P%^%s.prime)[current.state,2]
        }
    }
    return(var.total/window.length^2)
}


bar.sigmasq.fn <- function(P,pi, P.treat,window.length) {
    ## Estimate bar(sigma)^2
    eig.P.treat = eigen(P.treat)

    pi.treat = (eig.P.treat$vectors%*%diag(c(1,0))%*%solve(eig.P.treat$vectors))[1,]  # Stationary dist
    return(
    3/10*(pi.treat[1]*sigmasq.fn(P.treat,1,window.length) + pi.treat[2]*sigmasq.fn(P.treat,2,window.length)) +
    7/10*(pi[1]*sigmasq.fn(P,1,window.length) + pi[2]*sigmasq.fn(P,2,window.length))
    )
}
