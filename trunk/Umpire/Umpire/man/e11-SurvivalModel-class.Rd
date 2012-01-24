 \name{SurvivalModel-class}
\alias{SurvivalModel}
\alias{SurvivalModel-class}
\alias{rand,SurvivalModel-method}

\title{The "SurvivalModel" Class}
\description{A \code{SurvivalModel} object contains information for simulating survival times for each patient}
\usage{
SurvivalModel(baseHazard=1/5,
              accrual=5,
              followUp=1,
              units=12,
              unitName="months")
}
\arguments{
  \item{baseHazard}{numeric scalar describing the underlying hazard rate at baseline levels of covariates}
  \item{accrual}{numeric scalar representing number of patient accrual years}
  \item{followUp}{numeric scalar representing frequency of follow up in the unit of year}
  \item{units}{numeric scalar representing number of units per year where units are specified by \code{unitName}}
  \item{unitName}{character string representing the \code{unit} argument type}
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the
  \code{SurvivalModel} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{baseHazard}:}{Object of class \code{"numeric"}. see corresponding argument above}
    \item{\code{accrual}:}{Object of class \code{"numeric"}. see corresponding argument above}
    \item{\code{followUp}:}{Object of class \code{"numeric"}. see corresponding argument above}
    \item{\code{units}:}{Object of class \code{"numeric"}. see corresponding argument above}
    \item{\code{unitName}:}{Object of class \code{"character"}. see corresponding argument above}
  }
}
\section{Methods}{
  \describe{
    \item{rand(object, n, beta, \dots)}{simulate survival data for \code{n} patients given \code{beta}}
  }
}
\value{
  The \code{SurvivalModel} generator returns an object of class
  \code{SurvivalModel}.

  The \code{rand} method returns the survival data as a data.frame. The LFU column represents time to event, and Event column represents whether the data is censored.
}
\references{KRC}
\author{
  Kevin R. Coombes <kcoombes@mdanderson.org>
  Jiexin Zhang <jiexinzhang@mdanderson.org>
}
\seealso{
  \code{\linkS4class{CancerModel}}
}
\examples{
sm <- SurvivalModel()
## Generate data from base model
outcome <- rand(sm, 100)
summary(outcome)
## Generate data from five classes with different beta parameters
beta <- rep(rnorm(5, 0, 2), each=20)
outcome <- rand(sm, 100, beta=beta)
summary(outcome)
}
\keyword{classes}
\keyword{datagen}

