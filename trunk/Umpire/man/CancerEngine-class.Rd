\name{CancerEngine-class}
\docType{class}
\alias{CancerEngine-class}
\alias{rand,CancerEngine-method}
\alias{CancerEngine}

\title{The "CancerEngine" Class}
\description{A CancerEngine combines a CancerModel (which defines the
  combinatorics of hits that produce cancer subtypes) with a pair of
  gene expression Engines that can be used to simulate microarray data
  depending on the presence or absence of different hits.}
\usage{
CancerEngine(cm, base, altered)
}
\arguments{
    \item{cm}{A \code{\link{CancerModel}} object.}
    \item{base}{Object of class "character" giving the name of
      an \code{\link{Engine}} or \code{\link{EngineWithActivity}}, or an
      actual Engine.
      Represents the expected gene expression in the absence of any hits.}
    \item{altered}{Object of class "character" giving the name of
      an \code{\link{Engine}} or \code{\link{EngineWithActivity}},or an
      actual Engine.
      Represents the expected gene expression in the presence of hits.}
}
\section{Objects from the Class}{
Objects should be created using the function \code{CancerEngine}.
}
\section{Slots}{
	 \describe{
    \item{\code{cm}:}{A \code{\link{CancerModel}} object.}
    \item{\code{base}:}{Object of class "character" giving the name of
      an \code{\link{Engine}} or \code{\link{EngineWithActivity}}.
      Represents the expected gene expression in the absence of any hits.}
    \item{\code{altered}:}{Object of class "character" giving the name of
      an \code{\link{Engine}} or \code{\link{EngineWithActivity}}.
      Represents the expected gene expression in the presence of hits.}
    \item{\code{localenv}:}{Object of class \code{"environment"}; used
      in the internal implementation. }
  }
}
\section{Methods}{
  \describe{
    \item{rand}{\code{signature(object = "CancerEngine")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
\code{\linkS4class{CancerModel}}
}
\examples{
set.seed(391629)
# first build a CancerModel with 6 subtypes
nBlocks <- 20    # number of possible hits
cm <- CancerModel(name="cansim",
                  nPossible=nBlocks,
                  nPattern=6,
                  OUT = function(n) rnorm(n, 0, 1), 
                  SURV= function(n) rnorm(n, 0, 1))
# set up survival outcome; baseline is exponential
sm <- SurvivalModel(baseHazard=1/5, accrual=5, followUp=1)
# allow for activation of genes
activity <- rbinom(nBlocks, 1, 0.7)
# include blocks/pathways that are not hit by cancer
nExtraBlocks <- 100
nTotalBlocks <- nBlocks + nExtraBlocks
# assign values to hyperparameters
# block size
blockSize <- round(rnorm(nTotalBlocks, 100, 30))
# log normal mean hypers
mu0    <- 6
sigma0 <- 1.5
# log normal sigma hypers
rate   <- 28.11
shape  <- 44.25
# block corr
p <- 0.6
w <- 5
# set up the baseline Engine
rho <- rbeta(nTotalBlocks, p*w, (1-p)*w)
base <- lapply(1:nTotalBlocks, function(i) {
  bs <- blockSize[i]
  co <- matrix(rho[i], nrow=bs, ncol=bs)
  diag(co) <- 1
  mu <- rnorm(bs, mu0, sigma0)
  sigma <- matrix(1/rgamma(bs, rate=rate, shape=shape), nrow=1)
  covo <- co *(t(sigma) \%*\% sigma)
  MVN(mu, covo)
})
eng <- Engine(base)
# alter the means if there is a hit
altered <- alterMean(eng, normalOffset, delta=0, sigma=1)
# build the CancerEngine using character strings
object <- CancerEngine(cm, "eng", "altered")
# or build it using the actual Engine components
ob <- CancerEngine(cm, eng, altered)
# simulate some patients
clinical <- rand(cm, 20, sm)
dim(clinical)
summary(clinical)
# simulate some array data for these patients
class(object)
dset <- rand(object, clinical$CancerSubType)
summary(dset)
}
\keyword{classes}
\keyword{datagen}
