\name{CancerModel-class}
\docType{class}
\alias{CancerModel-class}
\alias{ncol,CancerModel-method}
\alias{nrow,CancerModel-method}
\alias{rand,CancerModel-method}
\alias{summary,CancerModel-method}

\title{Class "CancerModel" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
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
