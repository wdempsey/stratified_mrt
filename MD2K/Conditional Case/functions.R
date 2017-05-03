### Functions associated with Sample Size Calculations
require(foreach); require(TTR); require(zoo); require(expm)

### Simulation functions
rand.probs <- function(X.t, H.t, T, N, pi, lambda, min.p, max.p) {
  ## Calculate randomization probabilities given pi, x.t, lambda,
  ## T, and N
  power = length(H.t$X):1
  remaining.time = T-(max(power)+1)
  	if(remaining.time < 0) {
  		true.rem.time = 0
  	} else if(remaining.time - N[X.t]*60 < 60) {
  	  true.rem.time = remaining.time*pi[X.t]
  	} else if(remaining.time - N[X.t]*60 < 120) {
  		true.rem.time = (remaining.time-60)*pi[X.t]
	  } else {
		  true.rem.time = (remaining.time - 120)*pi[X.t]
	  }
	rho.t = max(min((N[X.t] - sum(((1-lambda^power)*H.t$rho + lambda^power*H.t$A)*(H.t$X == X.t & H.t$I == 1)))/ (1 + true.rem.time),max.p),min.p)
  return(rho.t)
}

calculateP <- function(inputs) {
    P.prime = matrix(0, nrow = 6, ncol = 6)
    P.prime[1,1] = inputs[1]; P.prime[4,4] = inputs[4]
    P.prime[2,3] = P.prime[5,6] = 1.0
    P.prime[1,2] = 1-P.prime[1,1]; P.prime[4,5] = 1-P.prime[4,4]
    P.prime[3,3] = inputs[2]; P.prime[6,6] = inputs[5]
    P.prime[3,1] = (1-inputs[3]) * (1-P.prime[3,3]); P.prime[3,4] = inputs[3] * (1-P.prime[3,3])
    P.prime[6,1] = (1-inputs[6]) * (1-P.prime[6,6]); P.prime[6,4] = inputs[6] * (1-P.prime[6,6])
    return(P.prime)
}

effect.gap <- function(P, window.length, effect) {
    f2 <- function(inputs) {
        P.prime = calculateP(inputs)
        total = 0.0
        for (k in 1:window.length) {
            total = total + (rowSums((P.prime%^%k)[c(2,5),4:6]) -
                rowSums((P%^%k)[c(2,5),4:6]))
        }
        gap = sum((total - effect*window.length)^2)
        return(gap)
    }
    return(f2)
}

test = optim(init.inputs,effect.gap(P, window.length, effect))
calculateP(test$par)

daily.sim <- function(N, pi, P.0, P.treat, T, window.length, min.p, max.p) {
  ## Simulate the Markov chain given length (T), transition matrix (P),
  ## and initial point (init)
  X.t = I.t = A.t = rho.t = vector(length = T)
  X.t[1] = sample(1:length(pi), size = 1, prob=pi)
  I.t[1] = as.numeric(X.t[1] == 2 | X.t[1] == 5)
  if(I.t[1] == 1) {
    rho.t[1] = max(min(N[X.t[1]]/((1 + (T - 120 - 1)*pi[X.t])),max.p),min.p)
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
            I.t[t] = as.numeric(X.t[t] == 2 | X.t[t] == 5)
            if( I.t[t] == 1) {
                rho.t[t] = rand.probs(X.t[t], H.t, T, N, pi, lambda, min.p, max.p)
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

sample.size <- function(ss.param,p,q,alpha.0 = 0.05,beta.0 = 0.8, max.iters = 10000) {
  ## Sample size calculation given the ss.param, p, and q

  p.calc <- function(N) {
    c.N = N*ss.param

    # inv.q  = (N-q-p)*(1-alpha.0)/ (p*(N-q-1))
    df1 = p
    df2 = N-q-p

    inv.f = (N - q - 1) * qf(1-alpha.0,df1,df2) / (N - q - p)

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

mean.Y <- function(P, window.length) {
    direct.result = 0
    for (k in 1:window.length) {
        direct.result = direct.result + (P%^%k)[,2]
    }

    return(direct.result/window.length)

}

ss.daily.data <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p, treatment.data){
  # Generate the daily data for a participant given all the inputs!

  inside.fn <- function(day) {
    effect = rep(daily.treat[day],2)
    P.treat = calc.Ptreat(P.0,effect,treatment.data, tol=10^(-2))
    H.t = daily.sim(N, pi, tau, P.0, P.treat, T, window.length, min.p, max.p)
    Y.t = SMA(H.t$X==2,window.length); Y.t = Y.t[(window.length+1):(length(Y.t))]
    prob.gamma = rollapply(1-H.t$rho, window.length, FUN = prod); prob.gamma = prob.gamma[-1]
    prob.nu = rollapply((H.t$A==0),window.length, FUN = prod); prob.nu = prob.nu[-1]
    psi.t = prob.nu/prob.gamma
    data = cbind(day,1:(T+window.length),H.t$A,H.t$X, H.t$rho,H.t$I, c(psi.t,rep(0,window.length)))
    Q = W = 0
    hat.sigmasq = list(); hat.sigmasq[[1]] = hat.sigmasq[[2]] = rep(0,0)
    hat.tildepr = hat.sigmasq
    obs.times = data[data[,6]==1 & data[,7] > 0,2]
    for(i in 1:(length(obs.times))) {
      ob.t = obs.times[i] ; f_at_obt = f.t(ob.t+ T*day, H.t$X[ob.t])
      p.t = tilde.p(H.t$X[ob.t])
      add.weight.t = (p.t/H.t$rho[ob.t])^(H.t$A[ob.t]) * ((1-p.t)/(1-H.t$rho[ob.t]))^(1-H.t$A[ob.t])
      P.ob.t = P*(1-H.t$A[ob.t]) + P.treat*H.t$A[ob.t]
      E.Y = mean.Y(P.ob.t, window.length)
      epsilon.ob.t = Y.t[ob.t] - E.Y[H.t$X[ob.t]]
      hat.sigmasq[[H.t$X[ob.t]]] = c(hat.sigmasq[[H.t$X[ob.t]]], (psi.t[ob.t] * add.weight.t * Y.t[ob.t] - E.Y[H.t$X[ob.t]])^2)
      hat.tildepr[[H.t$X[ob.t]]] = c(hat.tildepr[[H.t$X[ob.t]]], H.t$rho[ob.t])
      Q = Q +
          psi.t[ob.t] * add.weight.t * (H.t$A[ob.t]-p.t)^2 * outer(f_at_obt,f_at_obt)
#           p.t * (1-p.t) * outer(f_at_obt,f_at_obt)
      W = W +
          psi.t[ob.t] *
          add.weight.t *
          (H.t$A[ob.t]-p.t) *
          epsilon.ob.t *
          f_at_obt
    }
    return(list("Q.day" = Q, "W.day" = outer(W,W), "hat.sigmasq" = hat.sigmasq, "hat.tildepr" = hat.tildepr))
  }
  return(inside.fn)
}

full.trial.ss.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  W = Q = outer(f.t(1,1),f.t(1,1))*0
  hat.sigmasq = matrix(0,nrow = 2, ncol = 2)
  for(i in 1:length(daily.treat)) {
    day.res = ss.daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)(i)
    W = W + day.res$W.day
    Q = Q + day.res$Q.day
  }
  return(rbind(Q,W))
}

ss.parameters <- function(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  treatment.data = potential.effects(P, window.length)
  output = 0
  for (i in 1:num.iters) {
    output = output + full.trial.ss.sim(N, pi, tau,
                                        P.0, daily.treat, T,
                                        window.length, min.p, max.p,
                                        treatment.data)/num.iters
  }
  return(output)
}

####  Calculate bar.sigma

full.trial.barsigma.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  hat.sigmasq = matrix(0,nrow = 2, ncol = 2)
  for(i in 1:length(daily.treat)) {
    day.res = ss.daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)(i)
    hat.sigmasq = hat.sigmasq + rbind(c(sum(day.res$hat.sigmasq[[1]]), length(day.res$hat.sigmasq[[1]])), c(sum(day.res$hat.sigmasq[[2]]), length(day.res$hat.sigmasq[[2]])))
  }
  return(hat.sigmasq)
}

barsigma.estimation <- function(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  output = foreach(i=1:num.iters,.combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dopar% full.trial.barsigma.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)
  return(output)
}

# test.sigma = barsigma.estimation(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)


####  Calculate tilde.pr (optimal)

full.trial.tildepr.sim <- function(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data) {
  # Generate the full trial simulation using a vector of the daily treatment effects
  hat.tildepr = matrix(0,nrow = 2, ncol = 2)
  for(i in 1:length(daily.treat)) {
    day.res = ss.daily.data(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)(i)
    hat.tildepr = hat.tildepr + rbind(c(sum(day.res$hat.tildepr[[1]]), length(day.res$hat.tildepr[[1]])), c(sum(day.res$hat.tildepr[[2]]), length(day.res$hat.tildepr[[2]])))
  }
  return(hat.tildepr)
}

tildepr.estimation <- function(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p) {
  output = foreach(i=1:num.iters,.combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dopar% full.trial.tildepr.sim(N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p,treatment.data)
  return(output)
}

# test.tildepr = tildepr.estimation(num.iters, N, pi, tau, P.0, daily.treat, T, window.length, min.p, max.p)

## Sample Size simulation calculation

# Given Power Calc function and initial num.persons
# Compute the optimal number of people

power.calc <- function(num.persons, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data) {
  initial.study = foreach(k=1:1000, .combine = c,.packages = c('foreach','TTR','expm','zoo')) %dopar%
    estimation.simulation(num.persons, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
  return(mean(initial.study))
}

binary.search <- function(initial.N, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data) {

  High.N.current = High.N.old = round(initial.N*1.15,0)
  Mid.N.current = Mid.N.old = initial.N
  Low.N.current = Low.N.old = round(initial.N*0.85,0)

  which.run = rep(TRUE, 3); power.old = power.current = rep(0,3)
  total.evals = 0

  binary.iters = 100

  for(i in 1:binary.iters) {
    print(i)
    High.N.old = High.N.current
    Mid.N.old = Mid.N.current
    Low.N.old = Low.N.current
    power.old = power.current

    if(which.run[1] == TRUE) {
      power.L.current = power.calc(Low.N.current, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
    }
    if(which.run[2] == TRUE) {
      power.M.current = power.calc(Mid.N.current, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
    }
    if(which.run[3] == TRUE) {
      power.H.current = power.calc(High.N.current, N, pi, tau, P, daily.treat, T, window.length, min.p, max.p, treatment.data)
    }

    total.evals = total.evals + sum(which.run)

    if(power.L.current > power.M.current) {
      temp = mean(c(power.L.current,power.M.current))
      power.L.current = power.M.current = temp
    }

    if(power.M.current > power.H.current) {
      temp = mean(c(power.H.current,power.M.current))
      power.H.current = power.M.current = temp
    }

    if(power.L.current > power.H.current) {
      temp = mean(c(power.H.current,power.M.current,power.L.current))
      power.H.current = power.M.current = power.L.current = temp
    }

    if(power.L.current > 0.80) {
      Low.N.current = round(Low.N.current*0.8,0)
      Mid.N.current = round(Low.N.current*0.9,0)
      High.N.current = Low.N.current
      which.run = c(TRUE, TRUE, FALSE)
    } else if (power.H.current < 0.80) {
      High.N.current = round(High.N.current*1.2,0)
      Mid.N.current = round(High.N.current*1.1,0)
      Low.N.current = High.N.current
      which.run = c(FALSE, TRUE, TRUE)
    } else if (power.M.current < 0.80 ) {
      Mid.N.current = round(mean(c(Mid.N.current,High.N.current)),0)
      Low.N.current = Mid.N.old
      which.run = c(TRUE, FALSE, TRUE)
    } else {
      Mid.N.current = round(mean(c(Mid.N.current,Low.N.current)),0)
      High.N.current = Mid.N.old
      which.run = c(TRUE, FALSE, TRUE)
    }

    power.current = c(power.L.current,power.M.current,power.H.current)
    num.current = c(Low.N.current,Mid.N.current,High.N.current)

    print(rbind(num.current,power.current))
    print(total.evals)

    if(High.N.current - Mid.N.current <= 1 & Mid.N.current-Low.N.current <= 1) {
      best = min(which(power.current>0.80))
      return(c(num.current[best], power.current[best]))
    }

  }
  best = min(which(power.current>0.80))
  return(c(num.current[best], power.current[best]))
}
