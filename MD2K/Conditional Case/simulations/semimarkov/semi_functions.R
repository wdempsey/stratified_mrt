### Functions associated with Sample Size Calculations
require(foreach); require(TTR); require(zoo); require(expm); require(doRNG)

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

daily.sim <- function(N, pi, theta.0, theta.treat, T,
                      window.length, min.p, max.p,
                      pi.SMC) {
  ## Simulate the semi Markov chain given length (T),
  ## transition information
  states = state.list()
  I.t = A.t = rho.t = rep(0,length = T+2*window.length)
  X.t = U.t = L1.t = L2.t = L3.t = vector(length = T+2*window.length)
  init.state = states[sample(1:length(pi.SMC), size = 1, prob=pi.SMC),]
  t = 1
  X.t[t] = init.state[1]; U.t[t] = init.state[2]
  L1.t[t] = init.state[3]; L2.t[t] = init.state[4]; L3.t[t] = init.state[5]
  I.t[t] = as.numeric(U.t[t] == 2)
  if(I.t[t] == 1) {
    rho.t[t] = max(min(N[X.t[t]]/((1 + (T - 120 - 1)*pi.simple[X.t])),max.p),min.p)
  } else{rho.t[t] = 0}
  A.t[t] = rbinom(n=1,size=1, prob = rho.t[1])
  current.state = init.state
  while(t <= T+window.length) {
    if(A.t[t] == 1) {
      k = 0
      while(k  <=  window.length) {
        current.theta = theta.treat
        how.long = random.holding.time(current.state, current.theta)
        X.t[t:(t+how.long-1)] = current.state[1]
        U.t[t:(t+how.long-1)] = current.state[2]
        L1.t[t:(t+how.long-1)] = current.state[3]
        L2.t[t:(t+how.long-1)] = current.state[4]
        L3.t[t:(t+how.long-1)] = current.state[5]

        k = k + how.long
        t = t+how.long
        current.state = random.transition(current.state, current.theta)
      }
    } else{
      current.theta = theta.0
      how.long = random.holding.time(current.state, current.theta)
      X.t[t:(t+how.long-1)] = current.state[1]
      U.t[t:(t+how.long-1)] = current.state[2]
      L1.t[t:(t+how.long-1)] = current.state[3]
      L2.t[t:(t+how.long-1)] = current.state[4]
      L3.t[t:(t+how.long-1)] = current.state[5]
      current.state = random.transition(current.state, current.theta)
      t = t+how.long
    }
    X.t[t] = current.state[1]; U.t[t] = current.state[2]
    L1.t[t] = current.state[3]; L2.t[t] = current.state[4]
    L3.t[t] = current.state[5]
    I.t[t] = as.numeric(U.t[t] == 2) * (t <= T)
    H.t = list("X"=X.t[1:(t-1)],"A" = A.t[1:(t-1)],
               "I" = I.t[1:(t-1)], "rho" = rho.t[1:(t-1)])

    if(I.t[t] == 1) {
      rho.t[t] = rand.probs(X.t[t], H.t, T, N, pi, lambda, min.p, max.p)
    } else { rho.t[t] = 0 }
    A.t[t] = rbinom(n=1,size=1, prob = rho.t[t])
    if (any(is.na(H.t$A))) {print(t)}
  }
  return(H.t)
}

random.holding.time <- function(current.state, theta, max.hold = 100) {
  ## Generate random holding times from model given current state
  i = current.state
  if (i[2] == 2) {
    return(1)
  } else if (i[2] == 1) {
    est.shape = 1/theta$prepk.scale
    est.scale = exp(theta$prepk.coef%*%c(1,i[1]==2,i[3:5]==2))

    soj.prob = (pweibull(1:max.hold+1/2, scale = est.scale, shape = est.shape) -
                  pweibull(1:max.hold-1/2, scale = est.scale, shape = est.shape)) /
      (1-pweibull(1/2, scale = est.scale, shape = est.shape))

    return(
      sample(1:max.hold, size = 1, prob = soj.prob)
    )

  } else if (i[2] == 3) {
    est.shape = 1/theta$postpk.scale
    est.scale = exp(theta$postpk.coef%*%c(1,i[1]==2,i[3:5]==2))

    soj.prob = (pweibull(1:max.hold+1/2, scale = est.scale, shape = est.shape) -
                  pweibull(1:max.hold-1/2, scale = est.scale, shape = est.shape)) /
      (1-pweibull(1/2, scale = est.scale, shape = est.shape))

    return(
      sample(1:max.hold, size = 1, prob = soj.prob)
    )

  }
}

random.transition <- function(current.state, theta) {
  ## Generate random transition at from current state
  i = current.state
  states = state.list()
  comp.set = compatible.states(i)

  if (i[2] != 3) {
    return(
      as.numeric(comp.set)
    )
  } else {
    temp.trans = exp(theta$trans.coef%*%c(1,i[3:5]==2))
    prob.stress.trans = temp.trans/(1+temp.trans)

    temp.which = sample(1:2, size =1, prob = c(1-prob.stress.trans,
                                               prob.stress.trans))

    return(
      as.numeric(comp.set[,temp.which])
    )

  }

}

daily.data <- function(N, pi, theta.0, theta.treat,
                       T, window.length, min.p, max.p, pi.SMC){
  ## Generate the daily data for a participant given all the inputs.
  inside.fn <- function(day) {
    theta.treat = theta.treat.list[[day]]
    H.t = daily.sim(N, pi, theta.0, theta.treat, T,
                    window.length, min.p, max.p,
                    pi.SMC)
    Y.t = SMA(is.element(H.t$X,c(2)),window.length)
    Y.t = Y.t[(window.length+1):(T+window.length)]
    prob.gamma = rollapply(1-H.t$rho, window.length, FUN = prod)
    prob.gamma = prob.gamma[-1]
    prob.nu = rollapply((H.t$A==0),window.length, FUN = prod)
    prob.nu = prob.nu[-1]
    psi.t = prob.nu/prob.gamma
    data = cbind(day,1:T,Y.t,H.t$A[1:T],H.t$X[1:T],
                 H.t$rho[1:T],H.t$I[1:T], psi.t[1:T])
    return(data[data[,7] == 1 & data[,8] > 0,])
  }
  return(inside.fn)
}

full.trial.sim <- function(N, pi, theta.0, theta.treat.list,
                           T, window.length, min.p, max.p, pi.SMC) {
  ## Generate the full trial simulation for a person using a vector of the daily treatment effects
  foreach(i=1:length(theta.treat.list), .combine = "rbind",
          .packages = c("foreach", "TTR","expm","zoo")) %dorng% daily.data(N, pi, theta.0,
                                                                           theta.treat.list, T,
                                                                           window.length, min.p,
                                                                           max.p, pi.SMC)(i)
}

MRT.sim <- function(num.people, N, pi, theta.0, theta.treat.list,
                    T, window.length, min.p, max.p, pi.SMC) {
  ## Do the trial across people!!
  output = foreach(i=1:num.people, .combine = "rbind",
                   .packages = c("foreach", "TTR","expm","zoo")) %dorng% cbind(i,full.trial.sim(N, pi, theta.0, theta.treat.list,
                                                                                                T, window.length, min.p, max.p, pi.SMC))
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
  # Extract t(X)%*%W%*%X from the model
  X = cov[data[,1]==person,]
  W = exp(log.weights[data[,1]==person])
  t(X)%*%diag(W)%*%X
}


M.function <- function(cov, data, log.weights, person, XWX, fit) {
  # Compute the M term
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

  ## Set of possible weights depending on unique X.t
  set.rho = foreach(lvl=1:length(unique(X.t.person)), .combine = "c", .packages = c("foreach", "TTR","expm","zoo")) %dorng% mean(rho.t.person[X.t.person==lvl], na.rm = TRUE)

  rho = unlist(lapply(X.t.person,tilde.p))

  Z.t.person = B.t.person*matrix(rep(A.t.person-rho,3), ncol = 3)

  cov.t.person = cbind(B.t.person,Z.t.person)

  log.weights = A.t.person*(log(rho) - log(rho.t.person)) + (1-A.t.person)*(log(1-rho) - log(1-rho.t.person)) + log(psi.t.person)

  fit.people = lm(Y.t.person~(B.t.person+Z.t.person):as.factor(X.t.person)-1,weights = exp(log.weights))

  Covariates = model.matrix(fit.people)

  num.persons = length(unique(people[,1]))

  XWX = foreach(i=1:num.persons, .combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dorng% extract.tXWX(Covariates,people,log.weights,i)

  Middle = foreach(person=1:num.persons, .combine = "+", .packages = c("foreach", "TTR","expm","zoo")) %dorng% M.function(Covariates, people, log.weights, person,XWX,fit.people)

  entries1 = c(7:9)
  entries2 = c(10:12)
  entries = c(entries1,entries2)

  Sigma = solve(XWX,Middle)%*%solve(XWX)

  summ.fit.people = summary(fit.people)
  test.Sigma = summ.fit.people$sigma^2 * solve(crossprod(Covariates))


  output = (fit.people$coefficients[entries]%*%solve(Sigma[entries,entries],fit.people$coefficients[entries]))

  return(output)
}

estimation.simulation <- function(num.persons, N, pi, theta.0, theta.treat.list,
                                  T, window.length, min.p, max.p, pi.SMC) {

  people = MRT.sim(num.persons, N, pi, theta.0, theta.treat.list,
                   T, window.length, min.p, max.p, pi.SMC)

  output = estimation(people)

  alpha.0 = 0.05; p = 6; q = 6;

  multiple = p*(num.persons-q-1)/(num.persons-p-q)

  return (output>multiple*qf((1-alpha.0), df1 = p, df2 = num.persons - p - q))
}

#### SS-Calculation functions
tilde.p <- function(X.t) {
  ## Constant fn of t, and X.t
  ## N[2]/((T-60*N[2])*pi[2])*pi[2]+N[5]/((T-60*N[5])*pi[5])*pi[5]
  2*1.5/(600-60*1.5)
}

mean.Y <- function(P.ob.t, window.length) {
  total = 0.0
  for (k in 1:window.length) {
    total = total + (rowSums((P.ob.t%^%k)[c(2,5),4:6]))
  }
  return(c(0,total[1], 0, 0, total[2],0)/window.length)
}


### All functions related to the semiMarkov processes

compatible.states <- function (i) {
  ## Construct the set of states one
  ## can transition to from i
  if (i[2] == 1 | i[2] == 2) {
    j = as.numeric(i); j[2] = i[2] + 1
    return(as.matrix(j, ncol = 1))
  } else {
    j = as.numeric(i); j[2] = 1
    j[3] = i[1]; j[4] = i[3]; j[5] = i[4]
    return(as.matrix(cbind(c(1,j[2:5]), c(2,j[2:5])), nrow = 5))
  }
}

q_ij.k <- function(i, j, k, theta) {
  ## Function that returns probability of
  ## sojourn time = k given initially in
  ## state i = (X,U,L1,L2,L3)
  ## and then transitioning to state j

  comp.set = compatible.states(i)

  notcompatible = all(!(!colSums(comp.set != j ) ))

  if(notcompatible) {
    ## j is not compatible with i
    return( rep(0, length(k)) )
  } else if (i[2] == 1) {

    est.shape = 1/theta$prepk.scale
    est.scale = exp(theta$prepk.coef%*%c(1,i[1]==2,i[3:5]==2))

    soj.prob = (pweibull(k+1/2, scale = est.scale, shape = est.shape) -
                  pweibull(k-1/2, scale = est.scale, shape = est.shape)) /
      (1-pweibull(1/2, scale = est.scale, shape = est.shape))

    return(soj.prob)

  } else if (i[2] == 2) {

    return( as.numeric(k == 1) )

  } else if (i[2] == 3) {

    est.shape = 1/theta$postpk.scale
    est.scale = exp(theta$postpk.coef%*%c(1,i[1]==2,i[3:5]==2))


    soj.prob = (pweibull(k+1/2, scale = est.scale, shape = est.shape) -
                  pweibull(k-1/2, scale = est.scale, shape = est.shape)) /
      (1-pweibull(1/2, scale = est.scale, shape = est.shape))

    temp.trans = exp(theta$trans.coef%*%c(1,i[3:5]==2))
    prob.stress.trans = temp.trans/(1+temp.trans)

    prob.trans = prob.stress.trans*(j[1]==2) + (1-prob.stress.trans)*(j[1]==1)

    return(soj.prob * prob.trans)
  }

}

Q_ij.k <- function(i,j,k, theta) {
  ## Function that returns probability of
  ## sojourn time <= k given initially in
  ## state i = (X,U,L1,L2,L3)
  ## and then transitioning to state j

  sum(q_ij.k(i,j,1:k,theta))

}

Q_i.k <- function(i, k, theta) {

  comp.set = compatible.states(i)

  if(i[2] != 3) {
    j = as.numeric(comp.set)
    return ( Q_ij.k(i,j,k,theta) )
  } else{
    total = 0
    for (index in 1:2){
      j = as.numeric(comp.set[,index])
      total = total + Q_ij.k(i,j,k,theta)
    }
    return(total)
  }
}

h_i.k <- function(i,k, theta) {
  ## Function that returns the
  ## sojourn time distribution
  ## in state i

  comp.set = compatible.states(i)
  value = 0

  for (index in 1:ncol(comp.set)) {
    j = comp.set[,index]
    value = value + q_ij.k(i,j,k,theta)
  }

  return(value)

}

H_i.k <- function(i,k, theta) {
  ## Function that returns probability of
  ## sojourn time <= k given initially in
  ## state i = (X,U,L1,L2,L3)
  ## and then transitioning to state j

  sum(h_i.k(i,1:k,theta))

}

m_i <- function(i, theta, max.k = 100) {
  ## Function that returns
  ## the mean sojourn time in
  ## state i
  sum( (1:max.k)*h_i.k(i,1:max.k,theta) )
}

state.list <- function() {
  ## Construct the whole state list
  ## For the Markov chain
  a = c(1,2)
  b = c(1,2,3)

  as.matrix(expand.grid(a,b,a,a,a))
}

Markov.transitions <- function(theta) {
  ## Construct the Markov chain
  ## that ignores the holding times
  ## This is used to calculate the
  ## stationary distribution

  state = state.list()

  P = matrix(0, ncol = nrow(state), nrow = nrow(state))

  for(index in 1:nrow(P)) {

    i = state[index,]

    if(i[2] != 3) {
      j = as.numeric(compatible.states(i))

      which.column = which(!colSums(t(state) != j ) )

      P[index, which.column] = 1.0
    } else {
      comp.set = compatible.states(i)
      for (index.2 in 1:2) {
        j = as.numeric(comp.set[,index.2])
        which.column = which(!colSums(t(state) != j ) )

        temp.trans = exp(theta$trans.coef%*%c(1,i[3:5]==2))
        prob.stress.trans = temp.trans/(1+temp.trans)

        prob.trans = prob.stress.trans*(j[1]==2) + (1-prob.stress.trans)*(j[1]==1)
        P[index,which.column] = prob.trans
      }
    }
  }
  return(P)
}

stationary.dist <- function(theta) {
  P = Markov.transitions(theta)
  eig.P = eigen(P)

  pi = round(Re(eig.P$vectors%*%diag(c(1,rep(0,nrow(P)-1)))%*%solve(eig.P$vectors)),6)

  pi[1,]/sum(pi[1,])

}

limit.dist.SMC <- function(theta) {
  states = state.list()
  pi.SMC = rep(0, nrow(states))
  pi = stationary.dist(theta)
  for (index in 1:nrow(states)) {
    i = as.numeric(states[index,])
    pi.SMC[index] = m_i(i, theta) * pi[index]
  }
  pi.SMC/sum(pi.SMC)
}

p_ij.k <- function(i.loc, j.loc, k, states, theta, output) {
  ## Compute prob of being in state j
  ## after k steps from i
  ## Using past.array = old p_ij.k's

  i = as.numeric(states[i.loc,])
  j = as.numeric(states[j.loc,])

  if (k == 0) {
    return(as.numeric(all(i == j)))
  } else {

    comp.set = compatible.states(i)

    term1 = all(i==j) * (1- Q_i.k(i,k,theta))
    term2 = 0
    if (i[2] != 3) {
      l = as.numeric(comp.set)
      l.which.column = which(!colSums(t(states) != l ) )
      term2 = term2 + q_ij.k( i, l, 0:(k-1), theta) %*%
        output[l.which.column, j.loc, k:1]
    } else {
      for (index in 1:2) {
        l = as.numeric(comp.set[,index])
        l.which.column = which(!colSums(t(states) != l ) )
        term2 = term2 + q_ij.k( i, l, 0:(k-1), theta) %*%
          output[l.which.column, j.loc, k:1]
      }
    }
    return(term1 + term2)
  }

}

int.fn <- function(j.loc, k, states, theta, output) {
  i.s = 1:nrow(states)
  return(sapply(i.s, p_ij.k, j.loc, k, states, theta, output))
}

p_all.k <- function(Delta, theta) {
  ## Compute the probability
  ## of being in state j
  ## after k steps from
  ## step i for k = 1:Delta

  states <- state.list()
  num.states <- nrow(states)
  ## Output goes from step = 0 to Delta
  output = array(dim = c(num.states, num.states, Delta+1))
  for (k in 0:Delta) {
    output[,,k+1] = sapply(1:num.states, int.fn, k, states,theta,output)
  }

  return(output)
}

parallel.p_ij.k <- function(state.loc, k, states, theta, output) {
  ## Compute prob of being in state j
  ## after k steps from i
  ## Using past.array = old p_ij.k's

  i.loc = ceiling(state.loc/nrow(states))
  j.loc = state.loc - (i.loc - 1) * nrow(states)

  i = as.numeric(states[i.loc,])
  j = as.numeric(states[j.loc,])

  if (k == 0) {
    return(as.numeric(all(i == j)))
  } else {

    comp.set = compatible.states(i)

    term1 = all(i==j) * (1- Q_i.k(i,k,theta))
    term2 = 0
    if (i[2] != 3) {
      l = as.numeric(comp.set)
      l.which.column = which(!colSums(t(states) != l ) )
      term2 = term2 + q_ij.k( i, l, 0:(k-1), theta) %*%
        output[l.which.column, j.loc, k:1]
    } else {
      for (index in 1:2) {
        l = as.numeric(comp.set[,index])
        l.which.column = which(!colSums(t(states) != l ) )
        term2 = term2 + q_ij.k( i, l, 0:(k-1), theta) %*%
          output[l.which.column, j.loc, k:1]
      }
    }
    return(term1 + term2)
  }
}

parallel.p_all.k <- function(Delta, theta) {
  ## Compute the probability
  ## of being in state j
  ## after k steps from
  ## step i for k = 1:Delta

  states <- state.list()
  num.states <- nrow(states)
  ## Output goes from step = 0 to Delta
  output = array(dim = c(num.states, num.states, Delta+1))
  for (k in 0:Delta) {

    output[,,k+1] = matrix( parSapply(cl = cl, 1:num.states^2, parallel.p_ij.k, k=k, states=states,
                                      theta=theta, output = output),
                            nrow = num.states, ncol = num.states)

  }

  return(output)
}

proximal.outcome <- function(output, theta) {
  ## Returns expected fraction of time
  ## Stressed in next hour

  states <- state.list()

  stress.valid.states = as.logical(states[,1] == 2)

  fully.conditional.outcome = rowSums(output[,stress.valid.states,2:dim(output)[3]])/60

  pi.SMC = limit.dist.SMC(theta)

  return(c(
    sum(fully.conditional.outcome[!stress.valid.states]*
          pi.SMC[!stress.valid.states]/
          sum(pi.SMC[!stress.valid.states]))
    ,
    sum(fully.conditional.outcome[stress.valid.states]*
          pi.SMC[stress.valid.states]/
          sum(pi.SMC[stress.valid.states]))
  ))

}

treatment.effect<- function(baseline.prox, Delta,
                            alt.beta) {
  interior.fn <- function(theta.prime) {
    # print(theta.prime)

    theta.prime.df = list(
      "prepk.coef" = theta.prime[1:5],
      "prepk.scale" = theta.prime[6],
      "postpk.coef" = theta.prime[7:11],
      "postpk.scale" = theta.prime[12],
      "trans.coef" = theta.prime[13:16]
    )

    output.prime = p_all.k(Delta, theta.prime.df)
    treat.prox = proximal.outcome(output.prime, theta.prime.df)

    result = sum(
      60^2*((treat.prox - baseline.prox) - alt.beta)^2
    )

    return(result)

  }

  return(interior.fn)

}

optimal.treatment.day <- function(baseline.prox, Delta, daily.treat, day, init.theta) {
  if(abs(daily.treat[day]) < 10^-8) {
    day.barbeta.output = c(mean(daily.treat), day, unlist(temp.optim$par))

    return(init.theta)
  } else{
    alt.beta = rep(daily.treat[day], 2)

    treat.fn = treatment.effect(baseline.prox, Delta,
                                alt.beta)

    temp.optim = optim(init.theta,treat.fn, control = list(trace=TRUE, maxit = 2000))

    print(paste("Parameters =", temp.optim$par))
    print(paste("Value =", temp.optim$val))
    print(paste("Convergence =", temp.optim$convergence))

    print(paste("Done with day =", day, "for bar.beta =", round(mean(daily.treat),3)))

    temp.output = c(mean(daily.treat), day, unlist(temp.optim$par))

    day.barbeta.output = as.matrix(t(temp.output))

    write.table(day.barbeta.output, "export.csv", row.names = FALSE,
                col.names = FALSE, na = "NA", append = TRUE, sep = ",")

    temp.par = temp.optim$par
    saveRDS(temp.par, file = paste("output_",-mean(daily.treat),"_",day, ".rds", sep = ""))

    return(temp.optim$par)
  }
}

optimal.treatment.study <- function(baseline.prox, Delta, daily.treat, init.theta) {
  foreach(day=1:num.days, .combine = "cbind",
          .packages = c("foreach", "TTR","expm","zoo")) %dorng% optimal.treatment.day(baseline.prox, Delta, daily.treat, day, init.theta)
}

optimal.treatment.barbeta <- function(baseline.prox, Delta, bar.beta, init.theta) {
  ## Treatment vector
  Z.t = Vectorize(cov.gen)((1:num.days) * T)
  d = find.d(bar.beta,init.d,max.d,Z.t,num.days)
  daily.treat = -t(Z.t)%*%d

  return(optimal.treatment.study(baseline.prox, Delta, daily.treat, init.theta))

}

optimal.treatment.barbetaset <- function(baseline.prox, Delta, bar.beta.set, init.theta) {

  foreach(i=1:length(bar.beta.set), .combine = "cbind",
          .packages = c("foreach", "TTR","expm","zoo")) %dorng% optimal.treatment.barbeta(baseline.prox, Delta, bar.beta.set[i], init.theta)

}


relist.thetas <- function(unlisted.thetas) {
  ## Take the unlisted theta and make it into the
  ## list you need
  list.thetas = list(
    "prepk.coef" = unlisted.thetas[1:5],
    "prepk.scale" = unlisted.thetas[6],
    "postpk.coef" = unlisted.thetas[7:11],
    "postpk.scale" = unlisted.thetas[12],
    "trans.coef" = unlisted.thetas[13:16]
  )
  return(list.thetas)
}
