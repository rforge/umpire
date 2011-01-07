\name{blur-method}
\alias{blur}
\alias{blur-method}
\alias{blur,ANY-method}
\docType{methods}

\title{Methods blur}
\description{
  Functions to add noise to true signal.
}
\usage{
\S4method{blur}{NoiseModel}(object, x, \dots)
}

\arguments{
  \item{object}{object of class \code{NoiseModel}}
  \item{x}{The data matrix containing true signal from the gene expression}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\details{
  \code{blur} is the main method associated with a noise model.  The
    main operation is given by blur(object, x), which adds
  and multiplies random noise to the data matrix "x" containing the true signal
}
\references{KRC}
\author{Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}}
\seealso{
    NoiseModel
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
