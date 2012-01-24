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
  cancer patients.  It can be used to simulate patient outcome data
  linked to hit classes.}
\usage{
CancerModel(name, nPossible, nPattern,
                HIT = function(n) 5,
                SURV = function(n) rnorm(n, 0, 2),
                OUT = function(n) rnorm(n, 0, 2),
                survivalModel=NULL,
                prevalence=NULL)
nPatterns(object)
nPossibleHits(object)
nHitsPerPattern(object)
outcomeCoefficients(object)
survivalCoefficients(object)

}
\arguments{
  \item{name}{ The name given to this CancerModel. }
  \item{object}{ A CancerModel object. }
  \item{nPossible}{An integer; the number of potential hits relevant to
    the kind of cancer being modeled.}
  \item{nPattern}{The number of different cancer subtypes.}
  \item{HIT}{A function that generates non-negative integers from a
    discrete distribution.  Used to determine the number of hits
    actually present in each cancer subtype.}
  \item{SURV}{A function that generates real numbers from a continuos
    distributions.  Used in simulations to select the coefficients
    associated with each hit in Cox proportional hazards models.}
  \item{OUT}{A function that generates real numbers from a continuos
    distributions.  Used in simulations to select the coefficients
    associated with each hit in logistic models of a binary outcome.}
  \item{survivalModel}{An object of the SurvivalModel class.  used to
    simulate survival times for each simulated patient.}
  \item{prevalence}{ The prevalence of each cancer subtype. }
}

\section{Objects from the Class}{
Objects should be created using the \code{CancerModel} function.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"}}
    \item{\code{hitPattern}:}{Object of class \code{"matrix"}}
    \item{\code{survivalBeta}:}{Object of class \code{"numeric"}
      containing the coeffieicents associated with each hit in a
      Cox proportional hazards model of survival.}
    \item{\code{outcomeBeta}:}{Object of class \code{"numeric"}
      contaning the coefficients associated with each hit in a logistic
      model to predict a binary outcome.}
    \item{\code{prevalence}:}{Object of class \code{"numeric"}
      containing the prevalence of each cancer subtype. }
    \item{\code{survivalModel}:}{Object of class \code{"survivalModel"}
      containing parameters used to simualte survival times. }
    \item{\code{call}:}{Object of class \code{"call"} recording the
      function call used to initialize the object. }
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
  The \code{CancerModel} function is used to contruct and return an object of
  the \code{CancerModel} class.  The \code{ncol} and \code{nrow}
  functions return integers with the size of the matrix of hit
  patterns.  The \code{summary} function has no meanibngful return
  value; it is only called for its side-effects of printing a summary of
  the object.
  
  The \code{rand} method is the most important method for objects of this
  class.  It returns a data frame with four columns: the
  \code{CancerSubType} (as an integer that indexes into the
  \code{hitPattern} slot of the object), a binary \code{Outcome} that
  takes on values \code{"Bad"} or \code{"Good"}, an \code{LFU} column
  with censored survival times, and a logical \code{Event} column that
  describes whether the simulated survival event has occured.

  The \code{rand} method for the \code{CancerModel} class adds an extra
  logical parameter, \code{balanced}, to the signature specified by the
  default method.  If \code{balanced=FALSE} (the derfault), then
  patients are simulated based on the \code{prevalence} defined as apart
  of the model.   If \code{balanced=TRUE}, then patients are simulated
  with equal numbers in each hit pattern class, ordered by the hit
  pattern class.
}
\references{ OOMPA }
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
\code{\linkS4class{SurvivalModel}}
}
\examples{
set.seed(391629)
# set up survival outcome; baseline is exponential
sm <- SurvivalModel(baseHazard=1/5, accrual=5, followUp=1)
# now build a CancerModel with 6 subtypes
cm <- CancerModel(name="cansim",
                  nPossible=20,
                  nPattern=6,
                  OUT = function(n) rnorm(n, 0, 1), 
                  SURV= function(n) rnorm(n, 0, 1),
                  survivalModel=sm)
# simulate 100 patients
clinical <- rand(cm, 100)
summary(clinical)
}
\keyword{classes}
\keyword{datagen}
