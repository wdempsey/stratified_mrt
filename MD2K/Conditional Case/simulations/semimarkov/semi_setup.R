source("./semi_functions.R")

## Study assumptions

num.persons = 17 # Number of persons in the study

T = 600 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

## Duration coefficients
prepk.model <- readRDS("prepk_model.rds")
prepk.model$coefficients
prepk.model$scale

postpk.model <- readRDS("postpk_model.rds")
postpk.model$coefficients
postpk.model$scale

## Transition coefficients
transition.model <- readRDS("transition_model.rds")
transition.model$coefficients

## Aggregate the theta's into a list object
theta.0 <- list(
    "prepk.coef" = prepk.model$coefficients,
    "prepk.scale" = prepk.model$scale,
    "postpk.coef" = postpk.model$coefficients,
    "postpk.scale" = postpk.model$scale,
    "trans.coef" = transition.model$coefficients
)

## Construct the sequence of theta's
## in order to achieve the treatment effect
## on each day!

Delta = 60
outcome = p_all.k(Delta, theta)
baseline.prox = proximal.outcome(output)
theta.prime = unlist(theta) + rnorm(16, 0.001)
alt.beta = rep(-0.02,2)
temp.prime = treatment.effect(baseline.prox, Delta,
                              alt.beta)
treat.gap = temp.prime(theta.prime)
theta.vector = as.numeric(unlist(theta))
ptm <- proc.time()
test = optim(theta.vector, temp.prime, control = list(trace = TRUE))
proc.time() - ptm


# Randomization probability inputs
N =  c(0,1.61,0,0,2.05,0) # Avg. number of actions per day in each group ("Stressed", "Not Stressed")

lambda = 0.3 # Smoothing parameter in the randomization formula

min.p = 0.001 # Minimum randomization probability at each time point
max.p = 0.999 # Maximum randomization probability at each time point

## The stationary distribution for the simpler Markov example
## This is because it is what is used in the randomization
## algorithm.
pi.simple = c(0.0812, 0.0156)

# Treatment assumptions
init.d = 0 # Initial treatment effect
max.d = num.days/2 # Day of maximum treatment effect
bar.d = 0.01 # Avg treatment effect
