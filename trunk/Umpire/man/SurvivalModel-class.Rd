\name{SurvivalModel-class}
\alias{SurvivalModel}
\alias{SurvivalModel-class}

\title{Class "SurvivalModel" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SurvivalModel", ...)}.
	 ~~ describe objects here ~~ 
}
\usage{
SurvivalModel(baseHazard = 1/5, accrual = 5, followUp = 1, units = 12, unitName = "months")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseHazard}{ ~~Describe \code{baseHazard} here~~ }
  \item{accrual}{ ~~Describe \code{accrual} here~~ }
  \item{followUp}{ ~~Describe \code{followUp} here~~ }
  \item{units}{ ~~Describe \code{units} here~~ }
  \item{unitName}{ ~~Describe \code{unitName} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\section{Slots}{
	 \describe{
    \item{\code{baseHazard}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{accrual}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{followUp}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{units}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{unitName}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "SurvivalModel" in the signature.
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }
\seealso{~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{
showClass("SurvivalModel")
}
\keyword{classes}
\keyword{datagen}
