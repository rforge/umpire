\name{MVN-class}
\alias{MVN-class}
\alias{MVN}
\alias{alterMean,MVN-method}
\alias{alterSD,MVN-method}
\alias{nrow,MVN-method}
\alias{rand,MVN-method}
\alias{summary,MVN-method}
\alias{covar,MVN-method}
\alias{correl,MVN-method}

\title{The MV Class}
\description{The \code{MVN} class is a tool used to
  generate gene expressions that follow multivariate normal distribution}
\usage{
MVN(mu, Sigma, tol = 1e-06)
\S4method{alterMean}{MVN}(object, TRANSFORM, \dots)
\S4method{alterSD}{MVN}(object, TRANSFORM, \dots)
\S4method{nrow}{MVN}(x)
\S4method{rand}{MVN}(object, n, \dots)
\S4method{summary}{MVN}(object, \dots)
\S4method{covar}{MVN}(object)
\S4method{correl}{MVN}(object)
}
\arguments{
  \item{mu}{k-dimentional mean vector}
  \item{Sigma}{k*k covariance matrix containing the measurement of the
    linear coupling between every pair of randome vectors.}
  \item{tol}{Roundoff error that will be tolerated  when assessing the
    singularity of the covariance matrix}
  \item{object}{Object of class \code{MVN}}
  \item{x}{Object of class \code{MVN}}
  \item{TRANSFORM}{the \code{TRANSFORM} function for each object should
    take as its input a vector of mean expression or standard deviation
    and return a transformed vector that can be used to alter the
    appropriate slot of the object.}
  \item{n}{Number of samples to be simulated}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\section{Objects from the Class}{
Objects can be created by using the \code{MNV} generator function.  
}
\section{Slots}{
  \describe{
    \item{\code{mu}:}{Object of class \code{"numeric"} containing the k-dimentional mean vector}
    \item{\code{lambda}:}{Object of class \code{"numeric"} containing
      the square roots of eigenvalues of the covariance matrix}
    \item{\code{half}:}{a k*k matrix whose columns containing the
      eigenvectors of the covariance matrix}
  }
}

\section{Details}{
  The implementation of \code{MVN} class is designed for efficiency when
  generating new samples, since we expect to do this several times.
  Basically, this class separates the 'mvrnorm' function from the 'MASS'
  library into several steps.  The computationally expensive step (when
  the dimension is large) is the eigenvector decomposition of the
  covariance matrix. This step is performed at construction and the
  pieces are stored in the object. The \code{rand} method for \code{MVN} objects contains the second half
  of the 'mvrnorm' function from the 'MASS' library.
  
  Note that we typically work on expression value with its logarithm to
  some appropriate base. That is, the multivariate normal should
  be used on the logarithmic scale in order to contruct engine.

  \code{alterMean} for an \code{MVN} simply replaces the appropriate slot by
  the transformed vector. \code{alterSD} for an \code{MVN} is trickier, because of the way the data is
  stored. In order to have some hope of getting this correct, we
  work in the space of the covariance matrix, Sigma.  If we let
  R denote the correlation matrix and let Delta be the diagonal
  matrix whose entries are the individual standard deviations,
  then  Sigma = Delta %*% R %*% Delta.  So, we can change the
  standard deviations by replacing Delta in this product.  We then
  construct a new 'MVN' object with the old mean vector and the
  new covariance matrix.

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
    \item{alterMean(object, TRANSFORM,\dots)}{takes an object of class
      \code{MVN},loop over the \code{mu} slot, alter
      the mean as defined by TRANSFORM function, and returns an
      object of class \code{MVN} with altered \code{mu}}
    \item{alterSD(object, TRANSFORM,\dots)}{takes an object of class
      \code{MVN},works on the diagonal matrix of the covariance matrix, alter
      the standard deviation as defined by TRANSFORM function, and recontructs an
      object of class \code{MVN} with the old \code{mu} and
      reconstructed covariance matrix}
    \item{nrow(x)}{returns the number of genes (i.e, the
      length of the \code{mu} vector)}
    \item{rand(object, n, \dots)}{generates nrow(MVN)*n matrix representing gene
      expressions of n samples following the multivariate normal distribution
      captured in the object of \code{MVN}}
    \item{summary(object, \dots)}{prints out the number of
      multivariate normal random variables in the object of \code{MVN}}
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
    Engine, IndependentNormal
}
\examples{
\dontrun{
tolerance <- 1e-10
  # create a random orthogonal 2x2 matrix
  a <- runif(1)
  b <- sqrt(1-a^2)
  X <- matrix(c(a, b, -b, a), 2, 2)
  # now choos random positive squared-eigenvalues
  Lambda2 <- diag(rev(sort(rexp(2))), 2)
  # construct a covariance matrix
  Y <- t(X) %*% Lambda2 %*% X
  # Use the MVN constructor
  marvin <- MVN(c(0,0), Y)
  # check the four assertions
  print(paste('Tolerance for assertion checking:', tolerance))
  print(paste('Covar  assertion 1:',
              all(abs(covar(marvin) - Y) < tolerance)
              ))
  mar2 <- alterMean(marvin, normalOffset, delta=3)
  print(paste('Covar  assertion 2:',
              all(abs(covar(marvin) - covar(mar2)) < tolerance)
              ))
  print(paste('Correl assertion 1:',
              all(abs(diag(correl(marvin)) - 1) < tolerance)
              ))
  mar3 <- alterSD(marvin, function(x) 2*x)
  print(paste('Correl assertion 1:',
              all(abs(correl(marvin) - correl(mar2)) < tolerance)
              ))
  rm(a, b, X, Lambda2, Y, marvin, mar2, mar3)
}
}
\keyword{datagen}
\keyword{classes}
\keyword{distribution}
