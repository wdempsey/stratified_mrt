# Study assumptions

num.persons = 17 # Number of persons in the study

T = 600 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

## Construct Transition Matrix given Inputs 1 and 2
bar.W = c(0.074, 0.53)
bar.W.wkend1 = c(0.05, 0.50)
bar.W.wkend2 = c(0.10, 0.60)
bar.Z = c(10.3,11.8)

Pwkend.list = list(); piwkend.list = list()

source("./ti_functions.R")

P = alt.calculateP(bar.W,bar.Z)
pi = (P%^%1000)[1,]

P1 = alt.calculateP(bar.W.wkend1,bar.Z)
P2 = alt.calculateP(bar.W.wkend2,bar.Z)

Pwkend.list[[1]] = P1
Pwkend.list[[2]] = P2

approx.pi1 = (P1%^%1000)[1,]
approx.pi2 = (P2%^%1000)[1,]

piwkend.list[[1]] = approx.pi1
piwkend.list[[2]] = approx.pi2


# Randomization probability inputs
N =  c(0,1.61,0,0,2.05,0) # Avg. number of actions per day in each group ("Stressed", "Not Stressed")

lambda = 0.3 # Smoothing parameter in the randomization formula

min.p = 0.001 # Minimum randomization probability at each time point
max.p = 0.999 # Maximum randomization probability at each time point

# Treatment assumptions
init.d = 0 # Initial treatment effect
max.d = num.days/2 # Day of maximum treatment effect
bar.d = 0.01 # Avg treatment effect

## Initial inputs for finding the closest P.treat to give you the
## correct treatment effect.
init.inputs.wkend = list()
init.inputs1 = c(P1[1,1],P1[3,3],P1[3,4]/(1-P1[3,3]), P1[4,4], P1[6,6], P1[6,4]/(1-P1[6,6]))
init.inputs2 = c(P2[1,1],P2[3,3],P2[3,4]/(1-P2[3,3]), P2[4,4], P2[6,6], P2[6,4]/(1-P2[6,6]))

init.inputs.wkend[[1]] = init.inputs1
init.inputs.wkend[[2]] = init.inputs2

init.inputs.wkday = c(P[1,1],P[3,3],P[3,4]/(1-P[3,3]), P[4,4], P[6,6], P[6,4]/(1-P[6,6]))
