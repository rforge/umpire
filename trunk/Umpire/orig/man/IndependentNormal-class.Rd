\name{IndependentNormal-class}
\alias{IndependentNormal-class}
\alias{IndependentNormal}
\alias{alterMean,IndependentNormal-method}
\alias{alterSD,IndependentNormal-method}
\alias{nrow,IndependentNormal-method}
\alias{rand,IndependentNormal-method}
\alias{summary,IndependentNormal-method}

\title{The IndependentNormal Class}
\description{The \code{IndependentNormal} class is a tool used to
  generate gene expressions that follow independent normal distribution}
\usage{
IndependentNormal(mu,sigma)
\S4method{alterMean}{IndependentNormal}(object, TRANSFORM, \dots)
\S4method{alterSD}{IndependentNormal}(object, TRANSFORM, \dots)
\S4method{nrow}{IndependentNormal}(x)
\S4method{rand}{IndependentNormal}(object, n, \dots)
\S4method{summary}{IndependentNormal}(object, \dots)
}
\arguments{
  \item{mu}{Object of class \code{numeric} specifying the mean
    expression values.}
  \item{sigma}{Object of class \code{numeric} specifying the standard
    deviation of the gene expression values}
  \item{object, x}{Object of class \code{IndependentNormal}}
  \item{TRANSFORM}{the \code{TRANSFORM} function for each object should
    take as its input a vector of mean expression or standard deviation
    and return a transformed vector that can be used to alter the
    appropriate slot of the object.}
  \item{n}{Number of samples to be simulated}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\section{Objects from the Class}{
Objects can be created by using the \code{IndependentNormal} generator
function. The object of class \code{IndependentNormal} contains the mean
and standard deviation for the normal distribution
}
\section{Slots}{
	 \describe{
    \item{\code{mu}:}{Object of class \code{"numeric"} that contains the
    mean expression values}
    \item{\code{sigma}:}{Object of class \code{"numeric"} that contains
      the standard deviation of the gene expression values.}
  }
}
\section{Details}{
  Note that we typically work on expression value with its logarithm to
  some appropriate base. That is, the independent normal should
  be used on the logarithmic scale in order to contruct engine.
  }
\section{Methods}{
  \describe{
    \item{alterMean(object, TRANSFORM,\dots)}{takes an object of class
      \code{IndependentNormal},loop over the \code{mu} slot, alter
    the mean as defined by TRANSFORM function, and returns an
    object of class \code{IndependentNormal} with altered \code{mu}}
  \item{alterSD(object, TRANSFORM,\dots)}{takes an object of class
      \code{IndependentNormal},loop over the \code{sigma} slot, alter
    the standard deviation as defined by TRANSFORM function, and returns an
    object of class \code{IndependentNormal} with altered \code{sigma}}
    \item{nrow(x)}{\code{signature(object = "IndependentLogNormal")}: Returns the number of genes (i.e, the
      length of the \code{mu} vector)}
    \item{rand(object, n, \dots)}{generates nrow(IndependentNormal)*n matrix representing gene
      expressions of n samples following the normal distribution
      captured in the object of \code{IndependentNormal}}
    \item{summary(object, \dots)}{prints out the number of
      independent normal random variables in the object of \code{IndependentNormal}}
  }
}
\references{KRC}
\author{Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}}
\seealso{
  Engine, IndependentLogNormal, MVN
}
\examples{
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

}
\keyword{datagen}
\keyword{classes}
\keyword{distribution}
