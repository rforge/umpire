\name{CancerModel-class}
\alias{CancerModel-class}
\alias{CancerModel}
\alias{ncol,CancerModel-method}
\alias{nrow,CancerModel-method}
\alias{rand,CancerModel-method}
\alias{summary,CancerModel-method}
\alias{nPatterns}
\alias{nPossibleHits}
\alias{nHitsPerPattern}
\alias{survivalCoefficients}
\alias{outcomeCoefficients}

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
nPatterns(object)
nPossibleHits(object)
nHitsPerPattern(object)
outcomeCoefficients(object)
survivalCoefficients(object)

}
\arguments{
\item{name}{ ~~Describe \code{name} here~~ }
\item{object}{~~Describe \code{object} here~~ }
\item{NULL}{~~Describe \code{NULL} here~~ }
  \item{nPossible}{ ~~Describe \code{nPossible} here~~ }
  \item{nPattern}{ ~~Describe \code{nPattern} here~~ }
  \item{HIT}{ ~~Describe \code{HIT} here~~ }
  \item{SURV}{ ~~Describe \code{SURV} here~~ }
  \item{OUT}{ ~~Describe \code{OUT} here~~ }
  \item{prevalence}{ ~~Describe \code{prevalence} here~~ }
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
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }
\seealso{
\code{\linkS4class{SurvivalModel}}
}
\examples{
showClass("CancerModel")
}
\keyword{classes}
\keyword{datagen}
