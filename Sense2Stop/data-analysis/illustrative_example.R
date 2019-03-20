T = 6000
day_length = 600

day_in_study <- function(t) {
  floor((t-1)/day_length)+1
}

mu <- function(t) {
  beta_2 = -0.15
  beta_1 = (2-10)/10 - beta_2 * 10
  10 + beta_1 * (t/day_length) + beta_2 * (t/day_length)^2
}

Sigma <- function(t,s) {
  sigma_0 = 0.005
  sigma_1 = 0.5
  hour_corr = 0.8; hour_diff = 60
  lambda = (log(sigma_1) - log(sigma_1 + sigma_0) - log(hour_corr))/hour_diff
  sigma_0 * (t==s) + sigma_1 * exp( - lambda * abs(t-s))
}

Covariance_Matrix = matrix(nrow = T, ncol = T)

for(t in 1:T) {
  for (s in 1:T) {
    Covariance_Matrix[t,s] = Sigma(t,s)
  }
}

Chol = chol(Covariance_Matrix)

mean = vector(length = T)

for(t in 1:T) {
  mean[t] = mu(t)
}

set.seed(19783)
Z = rnorm(T)

outcome = t(Chol)%*%Z + mean

plot(1:T/day_length+1, outcome, type = "l")

# Construct Availability:
# Quadratic decrease from 0.8 to 0.2
exp_avail = vector(length = T)
for(t in 1:T) {
  # exp_avail[t] = 0.8-(0.8-0.2)/10 * t/day_length
  beta_2 = -0.01
  beta_1 = -beta_2 * 10 + (0.2-0.8)/10
  exp_avail[t] = 0.8 + beta_2 * (t/day_length)^2 + beta_1 *(t/day_length)
}

plot(exp_avail)

# Construct Marginal Prob of Stress:
# Quadratic from 0.15 up to 0.3 and back down to 0.15
exp_stress = vector(length = T)
for(t in 1:T) {
  # exp_avail[t] = 0.8-(0.8-0.2)/10 * t/day_length
  beta_2 = -0.006
  beta_1 = -beta_2 * 10 + (0.2-0.15)/10
  exp_stress[t] = 0.15 + beta_2 * (t/day_length)^2 + beta_1 *(t/day_length)
}

plot(exp_stress)


# Construct feature vectors
linear_f <- function(t) {
  return( c(1, day_in_study(t)) )
}

quadratic_f <- function(t) {
  return( c(1, day_in_study(t), day_in_study(t)^2 ) )
}


term1.quad = matrix(0, nrow = 3, ncol = 3)
term2.quad = rep(0,3)
term1.noavail.quad = matrix(0, nrow = 3, ncol = 3)
term2.noavail.quad = rep(0,3)
term1 = matrix(0, nrow = 2, ncol = 2)
term2 = rep(0,2)
term1.noavail = matrix(0, nrow = 2, ncol = 2)
term2.noavail = rep(0,2)
term1.const = matrix(0, nrow = 1, ncol = 1)
term2.const = rep(0,1)
term1.noavail.const = matrix(0, nrow = 1, ncol = 1)
term2.noavail.const = rep(0,1)

# Quad stress
term1.quad.quad = matrix(0, nrow = 3, ncol = 3)
term2.quad.quad = rep(0,3)
term1.linear.quad = matrix(0, nrow = 2, ncol = 2)
term2.linear.quad = rep(0,2)
term1.linear.noavail.quad = matrix(0, nrow = 2, ncol = 2)
term2.linear.noavail.quad = rep(0,2)
for(t in 1:T) {
  term2 = term2 + exp_avail[t]*outcome[t]* linear_f(t)
  term1 = term1 + exp_avail[t]*outer(linear_f(t),linear_f(t))
  term2.noavail = term2.noavail + outcome[t]* linear_f(t)
  term1.noavail = term1.noavail + outer(linear_f(t),linear_f(t))
  term2.const = term2.const + exp_avail[t]*outcome[t]
  term1.const = term1.const + exp_avail[t]
  term2.noavail.const = term2.noavail.const + outcome[t]
  term1.noavail.const = term1.noavail.const + 1
  term2.quad = term2.quad + exp_avail[t]*outcome[t]* quadratic_f(t)
  term1.quad = term1.quad + exp_avail[t]*outer(quadratic_f(t),quadratic_f(t))
  term2.noavail.quad = term2.noavail.quad + outcome[t]* quadratic_f(t)
  term1.noavail.quad = term1.noavail.quad + outer(quadratic_f(t),quadratic_f(t))
  
  term2.linear.quad = term2.linear.quad + exp_stress[t]*exp_avail[t]*outcome[t]* linear_f(t)
  term1.linear.quad = term1.linear.quad + exp_stress[t]*exp_avail[t]*outer(linear_f(t),linear_f(t))
  term2.linear.noavail.quad = term2.linear.noavail.quad + exp_stress[t]*outcome[t]* linear_f(t)
  term1.linear.noavail.quad = term1.linear.noavail.quad + exp_stress[t]*outer(linear_f(t),linear_f(t))
  term2.quad.quad = term2.quad.quad + exp_stress[t]*exp_avail[t]*outcome[t]* quadratic_f(t)
  term1.quad.quad = term1.quad.quad + exp_stress[t]*exp_avail[t]*outer(quadratic_f(t),quadratic_f(t))
}
beta.star = solve(term1,term2)
beta.star.noavail = solve(term1.noavail,term2.noavail)
beta.star.const = solve(term1.const,term2.const)
beta.star.noavail.const = solve(term1.noavail.const,term2.noavail.const)
beta.star.quad = solve(term1.quad,term2.quad)
beta.star.noavail.quad = solve(term1.noavail.quad,term2.noavail.quad)

beta.star.linear.quad = solve(term1.linear.quad,term2.linear.quad)
beta.star.linear.noavail.quad = solve(term1.linear.noavail.quad,term2.linear.noavail.quad)
beta.star.quad.quad = solve(term1.quad.quad,term2.quad.quad)




fitted_value = fitted_value.noavail = vector(length = T)
fitted_value.const = fitted_value.noavail.const = vector(length = T)
fitted_value.quad = fitted_value.noavail.quad = vector(length = T)

fitted_value.linear.quad = fitted_value.noavail.linear.quad = vector(length = T)
fitted_value.quad.quad = vector(length = T)

for(t in 1:T) {
  fitted_value[t] = linear_f(t)%*%beta.star
  fitted_value.noavail[t] = linear_f(t)%*%beta.star.noavail
  fitted_value.const[t] = beta.star.const
  fitted_value.noavail.const[t] = beta.star.noavail.const
  fitted_value.quad[t] = quadratic_f(t)%*%beta.star.quad
  fitted_value.noavail.quad[t] = quadratic_f(t)%*%beta.star.noavail.quad
  
  fitted_value.linear.quad[t] = linear_f(t)%*%beta.star.linear.quad
  fitted_value.noavail.linear.quad[t] = linear_f(t)%*%beta.star.linear.noavail.quad
  fitted_value.quad.quad[t] = quadratic_f(t)%*%beta.star.quad.quad

}

setwd("/Users/walterdempsey/Dropbox/MRTs/WD & SAM work/aoas MRTs/figs/")
png("intuition_figure.png", width = 480, height = 288, units = "px", pointsize = 12)
xseq = 1:T/day_length+1
par(mar = c(4.5,2,1,1)+0.1)
layout(matrix(c(1,1,1,1,1,1,1,1,2,2,3,3), nrow = 3, ncol = 4, byrow = TRUE))
plot(xseq, outcome, type = "l", col = "red", ylab="", xlab = "Day in Study",axes = FALSE)
axis(side=1, at = seq(1:11), labels = c(1:10,""))
axis(side=2, labels = FALSE)
points(xseq, fitted_value, type = "l", col = "blue")
points(xseq, fitted_value.noavail, type = "l", lty = 1, lwd = 3)
points(xseq, fitted_value.noavail.linear.quad, type = "l", col = "blue", lty=2)

points(xseq, fitted_value.noavail.quad, type = "l", lty = 3, lwd = 3)
#points(xseq, fitted_value.quad.quad, type = "l", lty = 3, col = "blue")

#points(xseq, fitted_value.const, type = "l", lty = 3)
#points(xseq, fitted_value.noavail.const, type = "l")

legend(1.5,7, c(expression(paste(beta, "(t;x)")), 
                expression(paste(f[t], "= (1,", d[t], "), quad E[", I[t], " | ", X[t],"= x], const P(",X[t], "=x)")), 
                expression(paste(f[t], "= (1,", d[t], "), const E[", I[t], " | ", X[t],"= x], const P(",X[t], "=x)")), 
                expression(paste(f[t], "= (1,", d[t], "), const E[", I[t], " | ", X[t],"= x], quad P(",X[t], "=x)")), 
                expression(paste(f[t], "= (1,", d[t], ", ", d[t]^2, "), const E[", I[t], " | ", X[t],"= x], const P(",X[t], "=x)"))),
       lty = c(1,1,1,2,3), lwd = c(1,3,1,1,3), 
       col = c("red", "black", "blue", "blue", "black"), bty = 'n')

plot(xseq, exp_avail, type = "l", col = "black", ylab="", ylim = c(0.0,1.0), xlab = "Day in Study",axes = FALSE)
lines(xseq, rep(mean(exp_avail), length(exp_avail)), col = "blue")
axis(side=1, at = seq(1:11), labels = c(1:10,""))
axis(side=2, at = seq(0.0,1.0,0.2), labels = c(0,rep("",4),1))
legend(1.5,0.7, c(expression(paste("Quadratic E[", I[t], " | ", X[t],"= x]")), 
                  expression(paste("Linear E[", I[t], " | ", X[t],"= x]"))),
       lty = c(1,1), lwd = c(1,1), 
       col = c("black", "blue"), bty = 'n')

plot(xseq, exp_stress, type = "l", col = "black", ylab="", ylim = c(0.0,1.0), xlab = "Day in Study",axes = FALSE)
lines(xseq, rep(mean(exp_stress), length(exp_stress)), col = "blue")
axis(side=1, at = seq(1:11), labels = c(1:10,""))
axis(side=2, at = seq(0.0,1.0,0.2), labels = c(0,rep("",4),1))
legend(1.5,1.1, c(expression(paste("Quadratic P(", X[t], "=x)")), 
                  expression(paste("Linear P(", X[t], "=x)"))),
       lty = c(1,1), lwd = c(1,1), 
       col = c("black", "blue"), bty = 'n')


dev.off()
