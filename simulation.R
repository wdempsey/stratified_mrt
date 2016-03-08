### Simulation setup
num.patients = 40 # number of participants
T = 12 * 60/5 # number decision points per day, 12 hours possible every five minutes
days = 14 # number of days
t = seq(1,days*T)
lambda = 1; N = c(2,2,0)

# For each t generate Z.t
cov.gen <-  function(t) {
  cov.t = c(1, floor((t-1)/T), floor((t-1)/T)^2)
  return(cov.t)
}

# Setup the treatment effect
bar.d = rep(avgtreatment,2)
init.d = rep(0,2)
max.d = 14/2
beta <- find.d(bar.d,init.d,max.d)
alpha <- rnorm(3,mean=0,sd=sd(beta)) 


## Generate action and response sequence

sim.patient <- function(rho) {
  X.t = A.t = rho.t = rep(0,0)
  for (i in 1:days) {
    init = sample(1:3, 1, replace = TRUE, prob= pi )
    daily.X.t = daily.sim(init,P,T)
    daily.rho.t = rand.probs(daily.X.t, T, N, pi, lambda, min.p, max.p)
    daily.A.t = d.rho.t = daily.X.t*0
    for (i in 1:T) {
      daily.A.t[i] = rbinom(1, size = 1, prob=daily.rho.t[daily.X.t[i],i])
      d.rho.t[i] = daily.rho.t[daily.X.t[i],i]
    }
    X.t = c(X.t, daily.X.t)
    A.t = c(A.t, daily.A.t)
    rho.t = c(rho.t, d.rho.t)
  }
  
  B.t = Z.t = Vectorize(cov.gen)(t)
  
  Y = t(B.t)%*%alpha + (A.t-rho)*t(Z.t)%*%beta + rnorm(length(t), 0, 1)
  
  return(list( "Y.t" = Y, "A.t" = A.t, "X.t" = X.t, "rho.t" = rho.t))
}

weight.fn <- function(A, rho, rho.hist) {
  res = exp(A* ( log(rho) - log(rho.hist) )  + (1-A)*(log(1-rho) - log(1-rho.hist)))
  res[is.nan(res)] = 0
  return(res)
}

update.est <- function(patient.data, rho) {
  W = weight.fn(patient.data$A.t, rho, patient.data$rho.t)
  B.t = Z.t = Vectorize(cov.gen)(t)
  Z.t = Z.t * (matrix(patient.data$A.t, nrow = nrow(Z.t), ncol = length(t), byrow=TRUE) - rho)  
  X = rbind(B.t, Z.t)
  
  num = den = list()
  for(k in 1:2) {
    num[[k]] = rep(0,nrow(X))
    den[[k]] = matrix(0, nrow = nrow(X), ncol = nrow(X))
  }
  
  for (i in t) {
    for(k in 1:2) {
      num[[k]] = num[[k]] + W[i]*patient.data$Y.t[i]*X[,i]*(patient.data$X.t[i]==k)
      den[[k]] = den[[k]] + W[i]*outer(X[,i],X[,i])*(patient.data$X.t[i]==k)
    }
  }
  return(list("num" = num, "den" = den))
}

beta.est <- function(num.patients, rho) {
  patient = num = den = list()
  for(k in 1:2) {
    num[[k]] = rep(0,nrow(X))
    den[[k]] = matrix(0, nrow = nrow(X), ncol = nrow(X))
  }
  
  for (pat in 1:num.patients) {
    patient[[pat]] <- sim.patient(rho)
    updates <- update.est(new.patient, rho)
    for (k in 1:2) {
      num[[k]] = num[[k]] + updates$num[[k]] 
      den[[k]] = den[[k]] + updates$den[[k]] 
    }  
  } 
  
  return(list("pat.data" = patient, "one" = solve(den[[1]],num[[1]]), "two" = solve(den[[2]],num[[2]])))
  
}

output = beta.est(num.patients,0.5)

cov.est <- function(output, rho) {
  
  for (pat in 1:num.patients) {
    rho.t = output$pat.data[[1]]$rho.t
    B.t = Z.t = Vectorize(cov.gen)(t)
    Z.t = Z.t * (matrix(output$pat.data[[1]]$A.t, nrow = nrow(Z.t), ncol = length(t), byrow=TRUE) - rho)  
    X = rbind(B.t,Z.t)
    for (i in 1:(T*days)) {
      if (rho > 0) {
        
        if (output$pat.data[[1]]$X.t[i] == 1) {
          beta = output$one
        } else if(output$pat.data[[1]]$X.t[i] ==2) {
          beta = output$two
        } else { beta= output$one * 0 }
      }
      output$pat.data[[1]]$Y.t[i] - X[,i]%*%beta*(output$pat.data[[1]]$A.t[i]-rho)
      
      
      B.t = Z.t = Vectorize(cov.gen)(t)
    }