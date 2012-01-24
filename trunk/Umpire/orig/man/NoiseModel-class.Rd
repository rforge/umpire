\name{NoiseModel-class}
\alias{NoiseModel-class}
\alias{NoiseModel}
\alias{blur,NoiseModel-method}

\title{The "NoiseModel" class}
\description{A \code{NOISE MODEL} represents the additional machine noise that is layered on top of any biological variabilty when measuring the gene expression in a set of samples}
\usage{
NoiseModel(nu, tau, phi)
\S4method{blur}{NoiseModel}(object, x, \dots)
}
\arguments{
  \item{nu}{The mean value for the additive noise}
  \item{tau}{The standard deviation for the additive noise}
  \item{phi}{The standard deviation for the multiplicative noise. Note
    the mean of multiplicative noise is set to 0}
  \item{object}{object of class \code{NoiseModel}}
  \item{x}{The data matrix containing true signal from the gene expression}
  \item{\dots}{extra arguments for generic or plotting routines}

}
\section{Details}{
  We model both additive and multiplicative noise, so that the observed
  expression of gene g in sample i is given by: Y_gi = S_gi exp(H_gi) +
  E_gi, where Y_gi = observed expression, S_gi = true bilogical signal,
  H_gi ~ N(0, phi) defines the multiplicative noise, and
  E_gi ~ N(nu,tau) defines the additive
  noise. Note that we allow a systematic offset/bias in the additive
  noise model.

  \code{blur} is the main method associated with a noise model.  The
  main operation is given by blur(object, x), which adds
  and multiplies random noise to the data matrix "x" containing the true signal.
}
\references{KRC}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}
}

\examples{

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
}
\keyword{classes}
\keyword{datagen}
