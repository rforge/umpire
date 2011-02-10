\name{CancerPatientSet-class}
\docType{class}
\alias{CancerPatientSet-class}
\alias{CancerPatientSet}
\alias{as.data.frame,CancerPatientSet-method}
\alias{summary,CancerPatientSet-method}

\title{Class "CancerPatientSet" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CancerPatientSet", ...)}.
	 ~~ describe objects here ~~ 
}
\usage{
CancerPatientSet(object, n)
\S4method{as.data.frame}{CancerPatientSet}(x, row.names = NULL, optional = FALSE)
}
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
    \item{row.names}{ ~~Describe \code{row.names} here~~ }
    \item{optional}{ ~~Describe \code{optional} here~~ }
}

\section{Slots}{
	 \describe{
    \item{\code{parent}:}{Object of class \code{"CancerModel"} ~~ }
    \item{\code{hitClass}:}{Object of class \code{"numeric"} ~~ }
  }
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "CancerPatientSet")}: ... }
    \item{summary}{\code{signature(object = "CancerPatientSet")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }
\seealso{
\code{\linkS4class{CancerModel}}
}
\examples{
showClass("CancerPatientSet")
}
\keyword{classes}
\keyword{datagen}
