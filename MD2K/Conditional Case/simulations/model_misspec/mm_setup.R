# Study assumptions

num.persons = 17 # Number of persons in the study

T = 600 # Number of timepoints per day

num.days = 10 # Number of days the study will be conducted

window.length = 60 # Window length for the calculating proximal outcome

# Construct Transition Matrix given Inputs 1 and 2
bar.W = c(0.074, 0.53)
bar.Z = c(10.3,11.8)

temp = c(1,-1)

epsilons = c(0.01,2)

potential.errors = expand.grid(temp,temp,temp,temp)

P.list = list(); pi.list = list()
temp.barW = bar.W
temp.barZ = bar.Z

for (i in 1:16) {
    temp.barW[1] = bar.W[1] + epsilons[1]*potential.errors[i,1]
    temp.barW[2] = bar.W[2] + epsilons[1]*potential.errors[i,2]

    temp.barZ[1] = bar.Z[1] + epsilons[2]*potential.errors[i,3]
    temp.barZ[2] = bar.Z[2] + epsilons[2]*potential.errors[i,4]

    P = alt.calculateP(temp.barW,temp.barZ)

    P.list[[i]] = P

    approx.pi = (P%^%1000)[1,]

    pi.list[[i]] = approx.pi

}

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
init.inputs = c(P[1,1],P[3,3],P[3,4]/(1-P[3,3]), P[4,4], P[6,6], P[6,4]/(1-P[6,6]))
