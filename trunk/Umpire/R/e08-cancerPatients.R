###############################################################
# The idea is that a CANCER MODEL can be used to generate a phenoData
# object for a simulated set of patients.  Combining this with a pair
# of ENGINE objects will later allow us to generate a complete exprSet
#
# Note that 'phenoData' and 'exprSet' are classes defined by BioConductor.

setClass("CancerPatientSet",
         representation=list(
           parent="CancerModel",
           hitClass="numeric"))

CancerPatientSet <- function(object, n) {
  if(!inherits(object, "CancerModel"))
    stop("First argument must be a 'CancerModel'")
  cp <- cumsum(object@prevalence)
  ru <- runif(n)
  hc <- unlist(lapply(ru, function(x, cp) {
    1 + length(cp) - sum(cp > x)
  }, cp)) # picks out a class by sampling from an explicit
          # discrete distribution
  new("CancerPatientSet",
      parent=object,
      hitClass=hc)
}

setMethod("summary", "CancerPatientSet", function(object, ...) {
  L <- length(object@hitClass)
  M <- max(object@hitClass)
  cat("A CancerPatientSet object representing outcome data\n")
  cat(paste("on", L, "patients from", M, "cancer subtypes,\n"))
  cat(paste('created from the cancer model: ',
            object@parent@name, '.\n', sep=''))
  X <- hist(object@hitClass, breaks=seq(0.5, M+0.5),
            plot=FALSE)$counts
  cat("Number of patients of each subtype:\n\n")
  print(X)
})

setMethod("as.data.frame", "CancerPatientSet",
          function(x, row.names=NULL, optional=FALSE) {
  data.frame()
})

realizeOutcome <- function(object) {
  if(!inherits(object, "CancerPatientSet"))
    stop("First argument must be a 'CancerPatientSet' object")
  hits <- object@parent
  temp <- as.vector(matrix(hits@outcomeBeta, nrow=1) %*%
                    hits@hitPattern)
  probs <- exp(temp)/(1+exp(temp))
  outclass <- c('Good', 'Bad')
  hc <- object@hitClass
  factor(outclass[1+rbinom(length(hc), 1, probs[hc])])
}

realizeSurvival <- function(object, baseHazard, accrual,
                            followUp=0, units=12) {
  if(!inherits(object, "CancerPatientSet"))
    stop("First argument must be a 'CancerPatientSet' object")
  hits <- object@parent
  temp <- as.vector(matrix(hits@survivalBeta, nrow=1) %*%
                    hits@hitPattern)
  hazard <- baseHazard*exp(temp)
  hc <- object@hitClass
  survival <- rexp(length(hc), hazard[hc]) # theoretical survival
  censor <- followUp + runif(length(hc), 0, accrual) # real time
  cen <- trunc(censor*units) # months, which is more likely ...
  sur <- trunc(survival*units)
  lfu <- ((sur + cen) - abs(sur-cen))/2 # vectorized minimum
  event <- (survival <= censor)         # did we observe the event?
  list(survival=lfu, event=event)
}

