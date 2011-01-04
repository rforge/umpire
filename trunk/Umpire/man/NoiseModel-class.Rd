\name{NoiseModel-class}
\alias{NoiseModel-class}
\alias{NoiseModel}
\alias{blur,NoiseModel-method}

\title{The "NoiseModel" class}
\description{A \code{NOISE MODEL} represents the additional machine noise that is layered on top of any biological variabilty when measuring the gene expression in a set of samples}
\usage{
NoiseModel(additiveOffset,additiveScale,multiplicativeScale)
\S4method{blur}{NoiseModel}(noiseModel, dataMatrix, \dots)
}
\arguments{
  \item{additiveOffset}{The mean value for the additive noise}
  \item{additiveScale}{The standard deviation for the additive noise}
  \item{multiplicativeScale}{The standard deviation for the multiplicative noise. Note
    the mean of multiplicative noise is set to 0}
  \item{noiseModel}{object of class \code{NoiseModel}}
  \item{dataMatrix}{The data matrix containing true signal from the gene expression}
}
\section{Details}
{
  We model both additive and multiplicative noise, so that the observed
  expression of gene g in sample i is given by: Y_gi = S_gi exp(H_gi) +
  E_gi, where Y_gi = observed expression, S_gi = true bilogical signal,
  H_gi ~ N(0, multiplicativeScale) defines the multiplicative noise, and
  E_gi ~ N(additiveOffset,additiveScale) defines the additive
  noise. Note that we allow a systematic offset/bias in the additive
  noise model.

  \code{blur} is the main method associated with a noise model.  The
  main operation is given by blur(moiseModel, dataMatrix), which adds
  and multiplies random noise to the dataMatrix containing the true signal.
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

additiveOffset <- 10
additiveScale <- 20
multiplicativeScale <- 0.1
nm <- NoiseModel(additiveOffset,additiveScale,multiplicativeScale)
realData <- blur(nm, myData)
summary(realData)
}
\keyword{classes}
\keyword{datagen}
