\name{covar-method}
\alias{correl}
\alias{correl, MVN-method}
\alias{covar}
\alias{covar, MVN-method}
\docType{methods}

\title{Methods covar and correl}
\description{
  Functions to extract covariance and correlation matrix specified in the \code{MVN} object.
}
\usage{
correl(object)
covar(object)
}

\arguments{
  \item{object}{object of class \code{MVN}}
}
\details{
  \code{covar} and \code{correl} functions calculate the covariance
    matrix and correlation matrix underlies the covariance matrix for the
    objects of class \code{MVN}, respectively. We
    have four assertions as shown below, and will be tested in the
    examples section:
    Assertion 1: \code{covar} should return the same matrix that was used
    in the function call to construct the MVN object.
    Assertion 2: After applying an "alterMean" function (below),
    the covariance matrix is unchanged.
    Assertion 3: The diagonal of correlation matrix consists of all 1's.
    Assertion 4: After applying an "alterMean" or an "alterSD"
  function (below), the correlation matrix is unchanged.
}
\section{Methods}{
  \describe{
    \item{covar(object)}{returns the covariance matrix of the object of
      class \code{MVN}}
    \item{correl(object)}{returns the correlation matrix of the object of
      class \code{MVN}}
  }
}
\references{KRC}
\author{Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}}
\seealso{
    MVN
}
