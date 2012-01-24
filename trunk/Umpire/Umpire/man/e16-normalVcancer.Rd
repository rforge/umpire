\name{NormalVsCancer}
\alias{NormalVsCancerModel}
\alias{NormalVsCancerEngine}

\title{Simulating Cancer Versus Normal Datasets}
\description{
  These functions are useful for simulating data that compares a
  homogeneous "cancer" group to a homogeneous "normal" group of
  samples.
}
\usage{
NormalVsCancerModel(nBlocks, name="NormalVsCancer")
NormalVsCancerEngine(nBlocks, hyperp) 
}
\arguments{
  \item{nBlocks}{An integer; the numbre of correlated blocks that are
    differentially expressed between cancer and normal.}
  \item{name}{A character string; the name of te object being created.}
  \item{hyperp}{An object of class \code{\link{BlockHyperParameters}}
    that describes the block correlation structure.}
}
\details{
  The simplest simulation model assumes that we are comparing two
  homogeneous groups.
}
\references{KRC}
\author{Kevin R. Coombes \email{kcoombes@mdanderson.org},
  Jiexin Zhang \email{jiexinzhang@mdanderson.org}}
\seealso{
  \code{\link{CancerModel}}, \code{\link{CancerEngine}},
  \code{\link{BlockHyperParameters}}.
}
\examples{
nvc <- NormalVsCancerModel(10)
summary(nvc)
rand(nvc, 10)
rand(nvc, 10, balance=TRUE)
engine <- NormalVsCancerEngine(10)
dset <- rand(engine, 10, balance=TRUE)
}
