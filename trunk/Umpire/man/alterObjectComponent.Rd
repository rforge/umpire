\name{alterObjectComponents-method}
\alias{alterMean}
\alias{alterSD}
\alias{normalOffset}
\alias{invGammaMultiple}
\docType{methods}
\title{Functions for alter components in the \code{Engine} object}
\description{
  \code{alterMean} and \code{alterSD} are functions that alter means or standard deviations in the
  \code{Engine} object. \code{normalOffset} is a possible
  \code{TRANSFORM} to be used in an \code{alterMean} operation, which
  adds an offset to each value in the mean where the offset is chosen
  from a normal distribution. \code{invGammaMultiple} is a possible
  \code{TRANSFORM} to be used in an \code{alterSD} operation, which
  multiplies each standard deviation by a positive value chosen from
  an inverse gamma distribution with paramters \code{shape} and \code{scale}.
}
\usage{
alterMean(object, TRANSFORM, ...)
alterSD(object,TRANSFORM,...)
normalOffset(x, delta, sigma)
invGammaMultiple(x, shape, rate)
}

\arguments{
  \item{object}{object of class \code{Engine}}
  \item{TRANSFORM}{the \code{TRANSFORM} function for each object should
    take as its input a vector of mean expression or standard deviation
    and return a transformed vector that can be used to alter the
    appropriate slot of the object.}
  \item{x}{numeric vector of mean expression or standard deviation
    defined in the \code{Engine} object}
  \item{delta, sigma}{mean and sd parameters specifying the normal distribution}
  \item{shape, rate}{shape and rate parameters specifiying the gamma distribution. Must be positive.}
  \item{\dots}{extra arguments for generic or plotting routines}
}

\value{
  \code{alterMean} and \code{alterSD} return a modified object of class
  \code{Engine} with corresponding components altered.

  \code{normalOffset} returns a new vector, each element of which is added by a offset
  chosen from a normal distribution with parameters \code{mean}
  and \code{sd}
  \code{invGammaMultiple} returns a new vector, each element of which is multiplied by a positive value
  chosen from an inverse gamma distribution with parameters \code{shape}
  and \code{scale}
  
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{distribution}
