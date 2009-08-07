\name{CancerModel}
\alias{CancerModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CancerModel(name, nPossible, nPattern, HIT = function(n) 5, SURV = function(n) rnorm(n, 0, 2), OUT = function(n) rnorm(n, 0, 2), prevalence = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{ ~~Describe \code{name} here~~ }
  \item{nPossible}{ ~~Describe \code{nPossible} here~~ }
  \item{nPattern}{ ~~Describe \code{nPattern} here~~ }
  \item{HIT}{ ~~Describe \code{HIT} here~~ }
  \item{SURV}{ ~~Describe \code{SURV} here~~ }
  \item{OUT}{ ~~Describe \code{OUT} here~~ }
  \item{prevalence}{ ~~Describe \code{prevalence} here~~ }
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
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
