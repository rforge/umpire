\name{SurvivalModel-class}
\alias{SurvivalModel}
\alias{SurvivalModel-class}
\alias{rand,SurvivalModel-method}

\title{The "SurvivalModel" Class}
\description{:TODO: class description}
\usage{
SurvivalModel(baseHazard=1/5,
              accrual=5,
              followUp=1,
              units=12,
              unitName="months")
}
\arguments{
  \item{baseHazard}{numeric scalar representing :TODO:}
  \item{accrual}{numeric scalar representing :TODO:}
  \item{followUp}{numeric scalar representing :TODO:}
  \item{units}{numeric scalar representing :TODO:}
  \item{unitName}{character string representing the \code{unit} argument type}
}
\details{
  :TODO: details
}
\section{Objects from the Class}{
  Although objects of the class can be created by a direct call to
  \link[methods]{new}, the preferred method is to use the
  \code{SurvivalModel} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{baseHazard}:}{see corresponding argument above}
    \item{\code{accrual}:}{see corresponding argument above}
    \item{\code{followUp}:}{see corresponding argument above}
    \item{\code{units}:}{see corresponding argument above}
    \item{\code{unitName}:}{see corresponding argument above}
  }
}
\section{Methods}{
  \describe{
    \item{rand(object, n, beta, \dots)}{:TODO:}
  }
}
\value{
  The \code{SurvivalModel} generator returns an object of class
  \code{SurvivalModel}.

  The \code{rand} method returns the survival data as a data.frame.
}
\references{KRC}
\author{
  Kevin R. Coombes <kcoombes@mdanderson.org>
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

