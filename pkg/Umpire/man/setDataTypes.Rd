\name{setDataTypes}
\alias{setDataTypes}
\title{
Discretize a Continuous Data Set to Mixed Types
}
\description{
The setDataTypes function allows the user to discretize a continuous data set 
generated from an engine of any type into binary, nominal, ordinal, or mixed-type data.
}
\usage{
setDataTypes(dataset, pCont, pBin, pCat, pNominal = 0.5, inputRowsAreFeatures = TRUE)
}
\arguments{
  \item{dataset}{A matrix or dataframe of continuous data.}
  \item{pCont}{
    Non-negative, numeric probability equal to or between 0 and 1 
    describing desired percent frequency of continuous features.
}
  \item{pBin}{
    Non-negative, numeric probability equal to or between 0 and 1 
    describing desired percent frequency of binary features.
}
  \item{pCat}{
    Non-negative, numeric probability equal to or between 0 and 1 
    describing desired percent frequency of categorical features.
}
  \item{pNominal}{
    Non-negative, numeric probability equal to or between 0 and 1 
    describing desired percent frequency of categorical features to 
    be simulated as nominal.
}
  \item{inputRowsAreFeatures}{
    Logical value indicating if features are to be simulated as rows 
    (TRUE) or columns (FALSE).
}
}
\details{:TO DO:
What level of detail should we go into here?
workflow; can store cutponts in mixed type engine
}
\value{
  \item{binned }{An object of class dataframe of simulated, mixed-type data.}
  \item{cutpoints }{
      An object of class "list". For each feature, data type, labels (if data are 
      binary or categorical), and break points (if data are binary or categorical) 
      are stored. Cutpoints can be passed to a \code{\link{MixedTypeEngine}} 
      for storage and downstream use.
      }
%% ...
}
\references{
Wang J, Wen S, Symmans WF, Pusztai L, Coombes KR. \cr\emph{The bimodality index: A criterion for discovering and ranking bimodal signatures from cancer gene expression profiling data.}\cr Cancer Inform 2009; 7:199-216.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
}
\note{
If \code{pCont}, \code{pBin}, and \code{pCat} do not sum to 1, they will be rescaled.

By default, categorical variables are simulates as a mixture of nominal and ordinal data. The remaining 
categorical features described by \code{pNominal} are simulated as ordinal.
}
\seealso{
  \code{\link{MixedTypeEngine}}
}
\examples{
## Generate a ClinicalEngine with 4 clusters
ce <- Umpire:::ClinicalEngine(20, 4, TRUE)

## Generate continuous data
set.seed(194718)
dset <- rand(ce, 300)

## Add noise before binning mixed type data
cnm <- ClinicalNoiseModel() # default
noisy <- blur(cnm, dset$data)

## Set the data mixture
dt <- setDataTypes(dset$data, 1/3, 1/3, 1/3, 0.3)
cp <- dt$cutpoints
type <- sapply(cp, function(X) { X$Type })
table(type)
summary(dt$binned)

## Use the pieces from above to create an MTE.
mte <- new("MixedTypeEngine",
           ce,
           noise = cnm,
           cutpoints = dt$cutpoints)
## and generate some data with the same data types and cutpoints
R <- rand(mte, 20)
summary(R)
}
\keyword{classes}
\keyword{datagen}
