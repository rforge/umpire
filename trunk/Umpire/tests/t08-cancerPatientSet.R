library(Umpire)
# set seed to ensure reprocuibility
set.seed(462662)
#  create a cancer model
cm <- CancerModel('simple',
                 30, 10,
                 HIT=function(n) 3+rbinom(1, 4, 0.7),
                 SURV=function(n) rnorm(n, 0, 1/2),
                 OUT=function(n) rnorm(n, 0, 2))
# generate a set of cancer patient from the cancer model
cps <- CancerPatientSet(cm, 100)
summary(cps)
# generate a set of outcomes for the patients
realizeOutcome(cps)
# generate survival data for the patients
temp <- data.frame(realizeSurvival(cps, baseHazard=1/5, accrual=5, units=52))
summary(temp)
