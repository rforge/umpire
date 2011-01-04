\name{CancerModel-class}
\alias{CancerModel-class}
\alias{CancerModel}
\alias{ncol,CancerModel-method}
\alias{nrow,CancerModel-method}
\alias{rand,CancerModel-method}
\alias{summary,CancerModel-method}

\title{The "CancerModel" Class}
\description{A \code{CancerModel} object contains a number of pieces of
  information representing an abstract, heterogeneous collection of
  cancer patients}
\usage{
CancerModel(name, nPossible, nPattern,
                HIT = function(n) 5,
                SURV = function(n) rnorm(n, 0, 2),
                OUT = function(n) rnorm(n, 0, 2),
                prevalence=NULL)
\S4method{alterMean}{Engine}(object, TRANSFORM, \dots)
\S4method{alterSD}{Engine}(object, TRANSFORM, \dots)
\S4method{nrow}{Engine}(x)
\S4method{rand}{Engine}(object, n, \dots)
\S4method{summary}{Engine}(object, \dots)
nComponents(object)
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CancerModel", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{hitPattern}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{survivalBeta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{outcomeBeta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{prevalence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{ncol}{\code{signature(x = "CancerModel")}: ... }
    \item{nrow}{\code{signature(x = "CancerModel")}: ... }
    \item{rand}{\code{signature(object = "CancerModel")}: ... }
    \item{summary}{\code{signature(object = "CancerModel")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }
\seealso{
}
\examples{
showClass("CancerModel")
}
\keyword{classes}
\keyword{datagen}
