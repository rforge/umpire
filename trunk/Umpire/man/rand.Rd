\name{rand-method}
\alias{rand}
\alias{rand-method}
\alias{rand,ANY-method}
\docType{methods}
\title{Method "rand"}
\description{
  Generate a matrix representing gene expressions following the
  distribution defined in the object.
}
\usage{
rand(object, n, ...)
}

\arguments{
  \item{object}{The \code{Engine} object defining the distributions of
    the gene expression}
  \item{n}{The number of samples being simulated}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\references{KRC}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}
}
\keyword{ datagen }
