library(lubridate)
library(plyr)

##### load data #####

dat <- read.csv("./data/stress_episodes.csv")
prob.dat <- read.csv("./data/stress_probabilities.csv")
lapse_timestamps <- read.csv("./data/lapse_timestamps.csv")
non_lapsers <- read.csv("./data/non_lapsers.csv")

### Timestamps for prob.dat ###

prob.dat$timestamp = as.POSIXct(prob.dat$Unix.Timestamp/1000, origin = "1970-01-01", "US/Central");

dat$st.time = as.POSIXct(dat$From.Timestamp.Unix./1000, origin = "1970-01-01", "US/Central");
dat$end.time = as.POSIXct(dat$To.Timestamp.Unix./1000, origin = "1970-01-01", "US/Central");
dat$peak.time = as.POSIXct(dat$Peak.Timestamp.Unix./1000, origin = "1970-01-01", "US/Central");


#### non_lapsers #####

lapse_timestamps <- subset(lapse_timestamps, Lapse.timestamp.Unix.!=-1) # remove the ones with no lapse timestamp
colnames(lapse_timestamps)[2] <- "relapseSession"
lapse_timestamps$relapse.time <- as.POSIXct(lapse_timestamps$Lapse.timestamp.Unix./1000, origin = "1970-01-01", "US/Central");
missing.data.id <- lapse_timestamps$Participant[!(lapse_timestamps$Participant %in% dat$Participant)]
lapse_timestamps <- subset(lapse_timestamps, !(Participant %in% missing.data.id))
non_lapsers$relapse <- 0

#### stress data cleaning####

wd.log <- dat;

# remove participants without relapse info
wd.log <- subset(wd.log, Participant %in% c(non_lapsers$Participant, lapse_timestamps$Participant))

# remove the episode with 0 duration in the begining of the day
wd.log <- subset(wd.log, end.time != st.time);

# Subset of data with start times before the peak time
wd.log <- subset(wd.log, st.time < peak.time)

# Compute durations and pre/post peak durations
wd.log$duration <- difftime(wd.log$end.time,wd.log$st.time, units = c("mins"));
wd.log$pre.pk.duration <- difftime(wd.log$peak.time,wd.log$st.time, units = c("mins"));
wd.log$post.pk.duration <- difftime(wd.log$end.time,wd.log$peak.time, units = c("mins"));

# Input 2 in Section 6.1.1
# For each episode type, the average episode lengt
mean(wd.log$duration[wd.log$stress == 1])  ## Episode = Not Stress
mean(wd.log$duration[wd.log$stress == 2])  ## Epsidoe = Stress

## Merge the two datasets
wd.log <- merge(wd.log, non_lapsers, by = c("Participant"), all=T)
wd.log <- merge(wd.log, lapse_timestamps, by=c("Participant"), all = T)

## remove pre-quit data
wd.log <- subset(wd.log, Session >10);
wd.log$Session <- wd.log$Session - 10;
wd.log$relapseSession <- wd.log$relapseSession - 10;

colnames(wd.log)[2] <- "day"
colnames(wd.log)[16] <- "relapse.day"

# no unsure class
stopifnot(all(wd.log$Class.1.NO.2.Unsure.3.YES.4.Unknown.!=2))

# rename dat columns
colnames(wd.log)[6] <- "stress.density"
colnames(wd.log)[7] <- "prop.missing"

# rename the stress
colnames(wd.log)[8] <- "stress"
wd.log$stress <- with(wd.log, (stress==1) + 2*(stress==3) + 3*(stress==4)) # Not Stress = 1, Stress = 2, Unsure = 3

# take the needed columns
wd.log <- subset(wd.log, select = c(Participant, day, stress.density, prop.missing, stress, st.time, end.time, duration, relapse.day, relapse, relapse.time, peak.time, pre.pk.duration, post.pk.duration))
wd.log$relapse <- with(wd.log, is.na(relapse) * 1)
wd.log$relapse.day <- with(wd.log, ifelse(is.na(relapse.day), 100, relapse.day))

# remove the insane unavailable
wd.log$st.hour <- as.numeric(format(wd.log$st.time, "%H"))
wd.log$end.hour <- as.numeric(format(wd.log$end.time, "%H"))

#stopifnot(log$stress[order(log$duration, decreasing = TRUE)[1:26]]==3)
index <- (wd.log$st.hour[order(wd.log$duration, decreasing = TRUE)[1:30]] < 6 | wd.log$st.hour[order(wd.log$duration, decreasing = TRUE)[1:30]] > 20)
wd.log <- wd.log[-order(wd.log$duration, decreasing = TRUE)[1:30][index], ]
wd.log <- wd.log[order(wd.log$Participant, wd.log$day, wd.log$st.time), ]


# redefine day by using 4am cutting point
define.day = function(log){
  log2 <- NULL
  for(id in unique(log$Participant)){
    
    dat <- subset(log, Participant == id)
    st.date <- as.POSIXct(format(dat$st.time[1], format="%Y/%m/%d"),  origin = "1970-01-01", "US/Central") + hours(4)
    index.date <- st.date + days(1) 
    
    temp.day <- as.numeric(with(dat, (end.time <= index.date)*1 + (end.time > index.date)* (1+ceiling(difftime(end.time, index.date, units = "hours")/24))))
    temp.day <- mapvalues(temp.day, sort(unique(temp.day)), 1:length(unique(temp.day)))
    
    dat$day <- temp.day
    
    #print(unique(dat$day))
    log2 <- rbind(log2, dat)
  }
  log2 <- log2[order(log2$Participant, log2$day, log2$st.time), ]
  return(log2)
}
wd.log <- define.day(wd.log)

# one user that never sleeps
wd.log <- subset(wd.log, Participant != 6032)

# redefine relapse day using the above cutting

define.lapseday = function(log){
  
  log2 <- NULL
  for(id in unique(log$Participant)){
    
    dat <- subset(log, Participant == id)
    
    if(is.na(dat$relapse.time[1])==F){
      dat$relapse.day <- dat$day[which.min(abs(as.numeric(difftime(dat$end.time, dat$relapse.time, units = "mins"))))]
    }
    
    log2 <- rbind(log2, dat)
  }
  return(log2)
}
wd.log = define.lapseday(wd.log)

### Discovering the classification rule from the stress density

# Proportion missing > 0.5 => Stress Classification "Unknown"

summary(wd.log$stress[wd.log$prop.missing > 0.5])
summary(wd.log$prop.missing[wd.log$stress == 3])

# Cutoff is 0.36 (below or equal = "No" and above = "Yes")
summary(wd.log$stress.density[(wd.log$prop.missing < 0.5) & (wd.log$stress == 1)])
summary(wd.log$stress.density[(wd.log$prop.missing < 0.5) & (wd.log$stress == 2)])


### Models for pre-pk, stress density, and then post-pk.  Each is an auto-regressive model where we use
### the past variables.  Based on the fitting of the models each ends up being maximum lag 2 or 1

# First we must construct robust summaries

require(foreign)
require(sandwich)

summary.robust <- function(model) {
  cov.m1 <- vcovHC(model, type = "HC0")
  
  std.err <- sqrt(diag(cov.m1))
  
  q.val <- qnorm(0.975)
  
  r.est <- cbind(
    Estimate = round(coef(model),2)
    , "Robust SE" = round(std.err,2)
    , z = round((coef(model)/std.err),2)
    , "Pr(>|z|) "= round(2 * pnorm(abs(coef(model)/std.err), lower.tail = FALSE),3)
    , LL = round(coef(model) - q.val  * std.err,2)
    , UL = round(coef(model) + q.val  * std.err,2)
  )
  
  r.est
}



auc.df = wd.log

library("dplyr")

auc.df <- 
  auc.df %>%
  group_by(Participant) %>%
  mutate(lag1.stress = dplyr::lag(stress, n = 1, default = NA),
         lag2.stress = dplyr::lag(stress, n = 2, default = NA),
         lag3.stress = dplyr::lag(stress, n = 3, default = NA),
         lag4.stress = dplyr::lag(stress, n = 4, default = NA),
         lag1.duration = dplyr::lag(duration, n = 1, default = NA),
         lag2.duration = dplyr::lag(duration, n = 2, default = NA),
         lag3.duration = dplyr::lag(duration, n = 3, default = NA),
         lag1.prepk = dplyr::lag(pre.pk.duration, n = 1, default = NA),
         lag2.prepk = dplyr::lag(pre.pk.duration, n = 2, default = NA),
         lag1.postpk = dplyr::lag(post.pk.duration, n = 1, default = NA),
         lag2.postpk = dplyr::lag(post.pk.duration, n = 2, default = NA))

## Input 1 in Section 6.1.1
ns.to.ns = sum(auc.df$stress[auc.df$lag1.stress == 1] == 1, na.rm = TRUE) # Number of observed transitions from "Not Stress" to "Not Stress"
ns.to.s = sum(auc.df$stress[auc.df$lag1.stress == 1] == 2, na.rm = TRUE) # Number of observed transitions from "Not Stress" to "Stress"

ns.to.ns/(ns.to.ns + ns.to.s)
ns.to.s/(ns.to.ns + ns.to.s)

s.to.ns = sum(auc.df$stress[auc.df$lag1.stress == 2] == 1, na.rm = TRUE) # Number of observed transitions from "Stress" to "Not Stress"
s.to.s = sum(auc.df$stress[auc.df$lag1.stress == 2] == 2, na.rm = TRUE) # Number of observed transitions from "Stress" to "Stress"

s.to.ns/(s.to.ns + s.to.s)
s.to.s/(s.to.ns + s.to.s)

### Stats on degree of missingness,
### What was the average number of hours of useable data per person?
users = auc.df$Participant
unique.users = unique(users)
diff.times = auc.df$end.time- auc.df$st.time
stress = auc.df$stress
lag1.stress = auc.df$lag1.stress

num.useable.hours = vector(length = length(unique.users))
num.useable.hours.lag1 = vector(length = length(unique.users))
num.useable.transitions.lag1 = vector(length = length(unique.users))
num.useable.stresstransitions = vector(length = length(unique.users))
for(id in 1:length(unique.users)) {
  num.useable.hours[id] = sum(as.numeric(diff.times[users == unique.users[id] & stress != 3]), na.rm = TRUE)/60
  num.useable.hours.lag1[id] = sum(as.numeric(diff.times[users == unique.users[id] & stress != 3 & lag1.stress != 3]), na.rm = TRUE)/60
  temp = as.numeric(diff.times[users == unique.users[id] & stress != 3 & lag1.stress != 3])
  num.useable.transitions.lag1[id] = length(temp[!is.na(temp)])
  temp = as.numeric(diff.times[users == unique.users[id] & (stress == 2 & lag1.stress == 2) ])
  num.useable.stresstransitions[id] = length(temp[!is.na(temp)])
}

num.useable.hours[is.na(num.useable.hours)] = 0.0
num.useable.hours.lag1[is.na(num.useable.hours.lag1)] = 0.0

hist(num.useable.hours)
mean(num.useable.hours)
sd(num.useable.hours)

hist(num.useable.hours.lag1)
mean(num.useable.hours.lag1)
sd(num.useable.hours.lag1)

mean(num.useable.transitions.lag1)
sd(num.useable.transitions.lag1)

mean(num.useable.stresstransitions)
sd(num.useable.stresstransitions)

### Fit the pre.pk model.  We check if the exponential distribution as a function
### of the stress classification is sufficient
library(survival)
auc.df$status = 1

baseline.model <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress), 
        data = subset(auc.df, (stress != 3) ), dist = "weibull", robust = TRUE)

summary(baseline.model)

model.lag1 <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress),
                      data = subset(auc.df, (stress != 3) & (lag1.stress != 3) ), dist = "weibull", robust = TRUE)

summary(model.lag1)

model.lag2 <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                        as.factor(lag2.stress) + as.numeric(lag1.postpk),
                      data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & (lag2.stress != 3) ), dist = "weibull", robust = TRUE)

summary(model.lag2)

model.lag3 <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                        as.factor(lag2.stress) + as.factor(lag3.stress),
                      data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & 
                                      (lag2.stress != 3) & (lag3.stress != 3) ), dist = "weibull", robust = TRUE)

summary(model.lag3)

complete.model.lag3 <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress)*as.factor(lag1.stress)*
                        as.factor(lag2.stress)*as.factor(lag3.stress),
                      data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & 
                                      (lag2.stress != 3) & (lag3.stress != 3) ), dist = "weibull", robust = TRUE)

test= summary(complete.model.lag3)$llik

2*(complete.model.lag3$loglik[2] - model.lag3$loglik[2]) > qchisq(0.95, df = 15 - 4)

model.lag4 <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                        as.factor(lag2.stress) + as.factor(lag3.stress) + as.factor(lag4.stress),
                      data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & 
                                      (lag2.stress != 3) & (lag3.stress != 3) & (lag4.stress != 3)), dist = "weibull", robust = TRUE)

summary(model.lag4)

# So the model is lag3 with lag1.postpk (length of prior window) # 
saveRDS(model.lag3,"prepk_model.rds")

## Fit the post.pk model.

postpk.baseline.model <- survreg(Surv(as.numeric(post.pk.duration), status) ~ as.factor(stress), 
                          data = subset(auc.df, (post.pk.duration > 0) & (stress != 3) ), dist = "weibull", robust = TRUE)

summary(postpk.baseline.model)

postpk.model.lag1 <- survreg(Surv(as.numeric(post.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress),
                      data = subset(auc.df, (post.pk.duration > 0) & (stress != 3) & (lag1.stress != 3) ), dist = "weibull", robust = TRUE)

summary(postpk.model.lag1)

postpk.model.lag2 <- survreg(Surv(as.numeric(post.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                               as.factor(lag2.stress),
                             data = subset(auc.df, (post.pk.duration > 0) & (stress != 3) & (lag1.stress != 3) & (lag2.stress != 3)), dist = "weibull", robust = TRUE)

summary(postpk.model.lag2)

postpk.model.lag3 <- survreg(Surv(as.numeric(post.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                               as.factor(lag2.stress)  + as.factor(lag3.stress),
                             data = subset(auc.df, (post.pk.duration > 0) & (stress != 3) & (lag1.stress != 3) 
                                           & (lag2.stress != 3) & (lag3.stress != 3)), dist = "weibull", robust = TRUE)

summary(postpk.model.lag3)

postpk.model.lag4 <- survreg(Surv(as.numeric(post.pk.duration), status) ~ as.factor(stress) + as.factor(lag1.stress) +
                               as.factor(lag2.stress) + as.factor(lag3.stress) + as.factor(lag4.stress),
                             data = subset(auc.df, (post.pk.duration > 0) & (stress != 3) & (lag1.stress != 3) 
                                           & (lag2.stress != 3) & (lag3.stress != 3) & (lag4.stress != 3)), dist = "weibull", robust = TRUE)

summary(postpk.model.lag4)


saveRDS(postpk.model.lag3,"postpk_model.rds")


### Now fit model for stress classification ###

stress.model.lag1 <- glm( I(stress == 2) ~ as.factor(lag1.stress), 
     data = subset(auc.df, (stress != 3) & (lag1.stress != 3) ),
     family = binomial)

summary.robust(stress.model.lag1)

stress.model.lag2 <- glm( I(stress == 2) ~ as.factor(lag1.stress) + as.factor(lag2.stress), 
                          data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & (lag2.stress != 3) ),
                          family = binomial)

summary.robust(stress.model.lag2)

stress.model.lag3 <- glm( I(stress == 2) ~ as.factor(lag1.stress) + 
                            as.factor(lag2.stress) + as.factor(lag3.stress), 
                          data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & (lag2.stress != 3) &
                                          (lag3.stress != 3)),
                          family = binomial)

summary.robust(stress.model.lag3)

saveRDS(stress.model.lag3,"transition_model.rds")

stress.model.lag4 <- glm( I(stress == 2) ~ as.factor(lag1.stress) + 
                            as.factor(lag2.stress) + as.factor(lag3.stress) + as.factor(lag4.stress), 
                          data = subset(auc.df, (stress != 3) & (lag1.stress != 3) & (lag2.stress != 3) &
                                          (lag3.stress != 3) & (lag4.stress != 3)),
                          family = binomial)

summary.robust(stress.model.lag4)


## Appears to be an AR(3) model ## 
## Here we calculate the transition rules ##
temp.stress.lag1 = exp(sum(stress.model.lag1$coefficients))
temp.notstress.lag1 = exp(stress.model.lag1$coefficients[1])

temp.stress.lag1/(1+temp.stress.lag1)
temp.notstress.lag1/(1+temp.notstress.lag1)

temp.allstress.lag3 = exp(sum(stress.model.lag3$coefficients))
temp.onlylag1stress.lag3 = exp(sum(stress.model.lag3$coefficients[1:2]))
temp.allnotstress.lag3 = exp(stress.model.lag3$coefficients[1])
temp.onlylag1notstress.lag3 = exp(sum(stress.model.lag3$coefficients[c(1,3:4)]))

temp.allstress.lag3/(1+temp.allstress.lag3)
temp.onlylag1stress.lag3/(1+temp.onlylag1stress.lag3)

temp.allnotstress.lag3/(1+temp.allnotstress.lag3)
temp.onlylag1notstress.lag3/(1+temp.onlylag1notstress.lag3)


## Marginal distribution plots for pre-peak and post-peak ##
library(ggplot2)

post.pk.df = subset(auc.df, (post.pk.duration > 0) & (post.pk.duration < 50) & (stress != 3) & (lag1.stress != 3) 
       & (lag2.stress != 3) & (lag3.stress != 3))

rep.num = length(as.numeric(post.pk.df$post.pk.duration))

post.pk.df$status = 1

## Frequency plots with fitted exponential and weibull curves
## 

weibull.plot.prepk.model <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ 1,
                          data = post.pk.df, dist = "weibull", robust = TRUE)

exp.plot.prepk.model <- survreg(Surv(as.numeric(pre.pk.duration), status) ~ 1,
                                    data = post.pk.df, dist = "exponential", robust = TRUE)

weibull.plot.postpk.model <- survreg(Surv(as.numeric(post.pk.duration), status) ~ 1,
                                    data = post.pk.df, dist = "weibull", robust = TRUE)

exp.plot.postpk.model <- survreg(Surv(as.numeric(post.pk.duration), status) ~ 1,
                                data = post.pk.df, dist = "exponential", robust = TRUE)

x.values = seq(0.5,40,0.5)
pre.y.exp.values = dweibull(x.values, shape = 1/exp.plot.prepk.model$scale, scale = exp(exp.plot.prepk.model$coefficients))
pre.y.wei.values = dweibull(x.values, shape = 1/weibull.plot.prepk.model$scale, scale = exp(weibull.plot.prepk.model$coefficients))

premodels.df = data.frame (x = x.values, y.exp = pre.y.exp.values, y.wei = pre.y.wei.values)


png("./figs/prepk_histogram.png", width = 480, height = 480, units = "px", pointsize = 12)
ggplot(post.pk.df, aes(as.numeric(pre.pk.duration))) + 
  geom_histogram(aes(y=..count../sum(..count..)), binwidth = 1, fill = "dodgerblue2", col = "dodgerblue") +
  geom_line(aes(x,y.exp), postmodels.df, col = "red", size = 1.25) +
  geom_line(aes(x,y.wei), postmodels.df, col = "black", size = 1.25) +
  ylab("Frequency") + 
  xlab("Duration (in minutes)")
dev.off()

x.values = seq(0.5,40,0.5)
post.y.exp.values = dweibull(x.values, shape = 1/exp.plot.postpk.model$scale, scale = exp(exp.plot.postpk.model$coefficients))
post.y.wei.values = dweibull(x.values, shape = 1/weibull.plot.prepk.model$scale, scale = exp(weibull.plot.postpk.model$coefficients))

postmodels.df = data.frame (x = x.values, y.exp = post.y.exp.values, y.wei = post.y.wei.values)

png("./figs/postpk_histogram.png", width = 480, height = 480, units = "px", pointsize = 12)
ggplot(post.pk.df, aes(as.numeric(post.pk.duration))) + 
  geom_histogram(aes(y=..count../sum(..count..)), binwidth = 1, fill = "dodgerblue2", col = "dodgerblue") +
  geom_line(aes(x,y.exp), postmodels.df, col = "red", size = 1.25) +
  geom_line(aes(x,y.wei), postmodels.df, col = "black", size = 1.25) +
  ylab("Frequency") + 
  xlab("Duration (in minutes)")
dev.off()
