pkgname <- "Umpire"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Umpire')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CancerEngine-class")
### * CancerEngine-class

flush(stderr()); flush(stdout())

### Name: CancerEngine-class
### Title: Class "CancerEngine" ~~~
### Aliases: CancerEngine-class rand,CancerEngine-method
### Keywords: classes datagen

### ** Examples

showClass("CancerEngine")



cleanEx()
nameEx("CancerModel-class")
### * CancerModel-class

flush(stderr()); flush(stdout())

### Name: CancerModel-class
### Title: The "CancerModel" Class
### Aliases: CancerModel-class CancerModel ncol,CancerModel-method
###   nrow,CancerModel-method rand,CancerModel-method
###   summary,CancerModel-method nPatterns nPossibleHits nHitsPerPattern
###   survivalCoefficients outcomeCoefficients
### Keywords: classes datagen

### ** Examples

showClass("CancerModel")



cleanEx()
nameEx("CancerPatientSet-class")
### * CancerPatientSet-class

flush(stderr()); flush(stdout())

### Name: CancerPatientSet-class
### Title: Class "CancerPatientSet" ~~~
### Aliases: CancerPatientSet-class CancerPatientSet
###   as.data.frame,CancerPatientSet-method summary,CancerPatientSet-method
### Keywords: classes datagen

### ** Examples

showClass("CancerPatientSet")



cleanEx()
nameEx("Engine-class")
### * Engine-class

flush(stderr()); flush(stdout())

### Name: Engine-class
### Title: The Engine Class
### Aliases: Engine-class Engine alterMean,Engine-method
###   alterSD,Engine-method nrow,Engine-method rand,Engine-method
###   summary,Engine-method nComponents
### Keywords: datagen classes

### ** Examples

nComp <- 10
nGenes <- 100
comp <- list()
for (i in 1:nComp){
  comp[[i]] <- IndependentNormal(rnorm(nGenes/nComp,6,1.5),1/rgamma(nGenes/nComp,44,28))
}
myEngine <- Engine(comp)
nrow(myEngine)
nComponents(myEngine)
summary(myEngine)
myData <- rand(myEngine,5)
dim(myData)
summary(myData)
OFFSET <- 2
myEngine.alterMean <- alterMean(myEngine,function(x){x+OFFSET})
myData.alterMean <- rand(myEngine.alterMean,5)
summary(myData.alterMean)
SCALE <- 2
myEngine.alterSD <- alterSD(myEngine,function(x){x*SCALE})
myData.alterSD <- rand(myEngine.alterSD,5)
summary(myData.alterSD)



cleanEx()
nameEx("EngineWithActivity-class")
### * EngineWithActivity-class

flush(stderr()); flush(stdout())

### Name: EngineWithActivity-class
### Title: The EngineWithActivity Class
### Aliases: EngineWithActivity-class EngineWithActivity
###   alterMean,EngineWithActivity-method alterSD,EngineWithActivity-method
###   rand,EngineWithActivity-method nrow,EngineWithActivity-method
###   summary,EngineWithActivity-method
### Keywords: datagen classes

### ** Examples

nComponents <- 10
nGenes <- 100
active <- 0.7
comp <- list()
for (i in 1:nComponents){
  comp[[i]] <- IndependentNormal(rnorm(nGenes/nComponents,6,1.5),1/rgamma(nGenes/nComponents,44,28))
}
myEngine <- EngineWithActivity(active,comp,2)
summary(myEngine)
myData <- rand(myEngine,5)
dim(myData)



cleanEx()
nameEx("IndependentLogNormal-class")
### * IndependentLogNormal-class

flush(stderr()); flush(stdout())

### Name: IndependentLogNormal-class
### Title: The IndependentNormal Class
### Aliases: IndependentLogNormal-class IndependentLogNormal
###   nrow,IndependentLogNormal-method rand,IndependentLogNormal-method
###   summary,IndependentLogNormal-method
### Keywords: datagen classes distribution

### ** Examples

  nGenes <- 20
  logmu <- rnorm(nGenes, 6, 1)
  logsigma <- 1/rgamma(nGenes, rate=14, shape=6)
  ln <- IndependentLogNormal(logmu, logsigma)
  nrow(ln)
  summary(ln)
  if(any(logmu - ln@logmu)) {
    print('means do not match')
  } else {
    print('means verified')
  }
  if(any(logsigma - ln@logsigma)) {
    print('standard deviations do not match')
  } else {
    print('sd verified')
  }
  x <- rand(ln, 1000)
  print(dim(x))
 
  print(paste("'ln' should be valid:", validObject(ln)))
  ln@logsigma <- 1:3 # now we break it
  print(paste("'ln' should not be valid:", validObject(ln, test=TRUE)))
  tmp.sd<-sqrt(apply(log(x),1,var))
  plot(tmp.sd,logsigma)
  tmp.mu<-apply(log(x),1,mean)
  plot(tmp.mu,logmu)
  rm(nGenes, logmu, logsigma, ln, x, tmp.mu, tmp.sd)



cleanEx()
nameEx("IndependentNormal-class")
### * IndependentNormal-class

flush(stderr()); flush(stdout())

### Name: IndependentNormal-class
### Title: The IndependentNormal Class
### Aliases: IndependentNormal-class IndependentNormal
###   alterMean,IndependentNormal-method alterSD,IndependentNormal-method
###   nrow,IndependentNormal-method rand,IndependentNormal-method
###   summary,IndependentNormal-method
### Keywords: datagen classes distribution

### ** Examples

nGenes <- 20
  mu <- rnorm(nGenes, 6, 1)
  sigma <- 1/rgamma(nGenes, rate=14, shape=6)
  ind <- IndependentNormal(mu, sigma)
  nrow(ind)
  summary(ind)
  if(any(mu - ind@mu)) {
    print('means do not match')
  } else {
    print('means verified')
  }
  if(any(sigma - ind@sigma)) {
    print('standard deviations do not match')
  } else {
    print('sd verified')
  }
  x <- rand(ind, 3)
  print(dim(x))
  print(summary(x))
  print(paste("'ind' should be valid:", validObject(ind)))
  ind@sigma <- 1:3 # now we break it
  print(paste("'ind' should not be valid:", validObject(ind, test=TRUE)))
  rm(nGenes, mu, sigma, ind, x)




cleanEx()
nameEx("MVN-class")
### * MVN-class

flush(stderr()); flush(stdout())

### Name: MVN-class
### Title: The MV Class
### Aliases: MVN-class MVN alterMean,MVN-method alterSD,MVN-method
###   nrow,MVN-method rand,MVN-method summary,MVN-method covar,MVN-method
###   correl,MVN-method
### Keywords: datagen classes distribution

### ** Examples

## Not run: 
##D tolerance <- 1e-10
##D   # create a random orthogonal 2x2 matrix
##D   a <- runif(1)
##D   b <- sqrt(1-a^2)
##D   X <- matrix(c(a, b, -b, a), 2, 2)
##D   # now choos random positive squared-eigenvalues
##D   Lambda2 <- diag(rev(sort(rexp(2))), 2)
##D   # construct a covariance matrix
##D   Y <- t(X) ##D 
##D   # Use the MVN constructor
##D   marvin <- MVN(c(0,0), Y)
##D   # check the four assertions
##D   print(paste('Tolerance for assertion checking:', tolerance))
##D   print(paste('Covar  assertion 1:',
##D               all(abs(covar(marvin) - Y) < tolerance)
##D               ))
##D   mar2 <- alterMean(marvin, normalOffset, delta=3)
##D   print(paste('Covar  assertion 2:',
##D               all(abs(covar(marvin) - covar(mar2)) < tolerance)
##D               ))
##D   print(paste('Correl assertion 1:',
##D               all(abs(diag(correl(marvin)) - 1) < tolerance)
##D               ))
##D   mar3 <- alterSD(marvin, function(x) 2*x)
##D   print(paste('Correl assertion 1:',
##D               all(abs(correl(marvin) - correl(mar2)) < tolerance)
##D               ))
##D   rm(a, b, X, Lambda2, Y, marvin, mar2, mar3)
## End(Not run)



cleanEx()
nameEx("NoiseModel-class")
### * NoiseModel-class

flush(stderr()); flush(stdout())

### Name: NoiseModel-class
### Title: The "NoiseModel" class
### Aliases: NoiseModel-class NoiseModel blur,NoiseModel-method
### Keywords: classes datagen

### ** Examples


nComp <- 10
nGenes <- 100
comp <- list()
for (i in 1:nComp){
  comp[[i]] <- IndependentLogNormal(rnorm(nGenes/nComp,6,1.5),1/rgamma(nGenes/nComp,44,28))
}
myEngine <- Engine(comp)
myData <- rand(myEngine,5)
summary(myData)

nu <- 10
tau <- 20
phi <- 0.1
nm <- NoiseModel(nu,tau,phi)
realData <- blur(nm, myData)
summary(realData)



cleanEx()
nameEx("SurvivalModel-class")
### * SurvivalModel-class

flush(stderr()); flush(stdout())

### Name: SurvivalModel-class
### Title: Class "SurvivalModel" ~~~
### Aliases: SurvivalModel SurvivalModel-class
### Keywords: classes datagen

### ** Examples

showClass("SurvivalModel")



cleanEx()
nameEx("alterObjectComponent")
### * alterObjectComponent

flush(stderr()); flush(stdout())

### Name: alterObjectComponents-method
### Title: Functions for alter components in the list("Engine") object
### Aliases: alterMean alterSD normalOffset invGammaMultiple
### Keywords: datagen distribution

### ** Examples

nComp <- 10
nGenes <- 100
comp <- list()
for (i in 1:nComp){
  comp[[i]] <- IndependentNormal(rnorm(nGenes/nComp,6,1.5),1/rgamma(nGenes/nComp,44,28))
}
myEngine <- Engine(comp)
nrow(myEngine)
nComponents(myEngine)
summary(myEngine)
myData <- rand(myEngine,5)
dim(myData)
summary(myData)
MEAN <- 2
SD <- 2
myEngine.alterMean <- alterMean(myEngine,function(x)normalOffset(x,
MEAN, SD))
myData.alterMean <- rand(myEngine.alterMean,5)
summary(myData.alterMean)
RATE <- 1
SHAPE <- 2
myEngine.alterSD <- alterSD(myEngine,function(x)invGammaMultiple(x,SHAPE,RATE))
myData.alterSD <- rand(myEngine.alterSD,5)
summary(myData.alterSD)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
