## Copyright (C) Kevin R. Coombes and Caitlin E. Coombes, 2020

## Objects of this class remember the entire procedure to
## (a) simulate raw (continuous) data,
## (b) add a particular kind of noise to it, and
## (c) discretize some features.
setClass("MixedTypeEngine",
         contains = "CancerEngine",
         slots = c(noise = "NoiseModel",
                   cutpoints = "list"))

## Given the desired percentages or fractions of each data type
## and a "training" data set, this function decides which features
## to make binary or categorical. It also records the cutpoints,
## associated labels, and type of each feature. The result can be
## used to create a "MixedTypeEngine" both for the ability to
## report parameters and to generate later test sets fromthe same
## multivariate distribution.
setDataTypes <- function(dataset, pCont, pBin, pCat, pNominal=0.5) {
  nFeats <- nrow(dataset)
  cutpoints <- list() # prepare to fill a list
  probs <- c(pCont, pBin, pCat)
  if (any(probs < 0)) stop("All type probabilities must be nonnegative.")
  if (sum(probs) != 1) {
    warning("Probabilities are being rescaled so they add to 1.")
    probs <- probs/sum(probs)
  }
  Type <- sample(c("continuous", "binary", "categorical"),
                 size = nFeats, replace = TRUE, prob = probs)
  isCont <- Type == "continuous"
  for (W in which(isCont)) {
    L <- list(breaks = NULL, labels = NULL, Type = "continuous")
    cutpoints[[W]] <- L
  }
  binned <- dataset
  if (any(isbin <- Type == "binary")) {
    X <- makeBinary(dataset[isbin, , drop = FALSE])
    binned[isbin,] <- X$dataset
    cutpoints[isbin] <- X$cuts
  }
  if (any(iscat <- Type == "categorical")) {
    X <- makeCategories(dataset[iscat, , drop = FALSE], pNominal)
    binned[iscat,] <- X$dataset
    cutpoints[iscat] <- X$cuts
  }
  list(binned = binned, cutpoints=cutpoints)
}


makeCategories <- function(dataset, pNominal = 0.5) {
  nPts <- ncol(dataset)
  cutpoints <- list() # prepare to fill a list
  ## loop through rows of dataset
  for (I in 1:nrow(dataset)) {
    ## decide if variable is ordinal or nominal
    catType <- sample(c("ordinal", "nominal"), 1, prob = c(1 - pNominal, pNominal))
    ## Sample a number of categories between 3 and 9.
    cc <- sample(3:9, 1, replace = TRUE)
    ## Create a numbered list of bin identities.
    if (catType == "ordinal") {
      id <- 1:cc
    } else { # nominal
      id <- sample(1:cc, cc, replace = FALSE)
    }
    ## Sample bin sizes from the Dirichlet distribution
    r <- rdirichlet(1, rep(20, cc)) # magic=20 determined empirically
    ## Construct a vector of percentage cutoffs
    pcuts <- c(0, cumsum(r))
    qcuts <- quantile(dataset[I,], pcuts)
    ## move the extremes out to avoid problems on future data
    qcuts[1] <- -Inf
    qcuts[length(qcuts)] <- Inf
    ## bin everything
    M <- cut(dataset[I,], breaks=qcuts, labels=id, include.lowest=TRUE)
    cutpoints[[I]] <- list(breaks  = qcuts,
                           labels = id,
                           Type = catType)
    dataset[I,] <- M
  } # end loop through rows of dataet
  list(dataset = dataset, cuts = cutpoints)
} # end function

makeBinary <- function(dataset) {
  ## KRC: May need to check that BI works on single-row data matrix
  ## if it looks bimodal, split at the midpoint of modes
  bi <- bimodalIndex(dataset) # use the bimodal index (Wang et al)
  my.bi <- 1.1                # and cutoff from the paper
  cutpoints <- list()
  for (I in 1:nrow(dataset)) {
    if(bi[I, "BI"] > my.bi){
      bound <- (bi[I,"mu1"] + bi[I,"mu2"])/2
    } else {
      perc <- sample(5:35, 1) # arbitrary split point
      perc <- perc/100       # converted to fraction
      # decide whether to keep lower or upper chunk
      low.cutoff <- quantile(dataset[I,], perc)
      high.cutoff <- quantile(dataset[I,], 1 - perc)
      bound <-  ifelse(bi[I,"pi"] < 0.5, low.cutoff, high.cutoff)
    }
    bks <- c(-Inf, bound, Inf)
    lbls <- sample(c(0,1)) # randomly assign 0, 1 values.
    X <- cut(dataset[I,], breaks = bks, labels = lbls)
    X <- as.numeric(as.character(X))
    type <- ifelse(mean(X) < 0.10, "asymmetric binary", "symmetric binary")
    cutpoints[[I]] <- list(breaks = bks, labels = lbls, Type = type)
    dataset[I,] <- X
  }
  list(dataset = dataset, cuts = cutpoints)
}


## It's probably worth some time to figure out all the properties
## that an MTE should really have, and put all those checks here..
setValidity("MixedTypeEngine", function(object) {
  msg <- NULL
  rightSize <- length(object@cutpoints) == nrow(object)
  cat("|", rightSize, "|, L =", length(rightSize), "\n", file=stderr())
#  if (!rightSize) {
#    msg <- c(msg, "Sizes of cut points and engine do not match.")
#  }
  ## make sure cutpoints are sensible
  listCheck <- all (sapply(object@cutpoints, is.list))
  if (!listCheck) {
    msg <- c(msg, "All elements of 'cutpoints' must themselves be lists.")
  } else {
    hasBreaks <- all (sapply(object@cutpoints, function(X) {
      "breaks" %in% names(X)
    }))
    hasLabels <- all (sapply(object@cutpoints, function(X) {
      "labels" %in% names(X)
    }))
    hasType <- all (sapply(object@cutpoints, function(X) {
      "Type" %in% names(X)
    }))
    if (!hasBreaks) {
      msg <- c(msg, "Every element of 'cutpoints' must include 'breaks'.")
    }
    if (!hasLabels) {
      msg <- c(msg, "Every element of 'cutpoints' must include 'labels'.")
    }
    if (!hasType) {
      msg <- c(msg, "Every element of 'cutpoints' must include the 'Type' of the feature.")
    }
    if (hasBreaks & hasLabels & hasType) {
      nullBreaks <- sapply(object@cutpoints, function(OC) is.null(OC$breaks))
      nullLabels <- sapply(object@cutpoints, function(OC) is.null(OC$labels))
      isCont <- sapply(object@cutpoints, function(OC) OC$Type == "continuous")
      cat(table(isCont), "\n", file = stderr())
      cat(table(nullBreaks), "\n", file = stderr())
      if (!all(isCont == nullBreaks)) {
        msg <- c(msg, "Null break list does not agree with 'continuous' type value.")
      }
      if (!all(nullBreaks == nullLabels)) {
        msg <- c(msg, "Null breaks does not match null labels.")
      }
    }
  }
  if (is.null(msg)) msg <- TRUE
  msg
})

## 'rand' generates random-vectors from a given joint distibution. This includes
## adding noise and discretizing consistently. Mostly used for test data, since
## the construction of a realistic MTE requires generating a training data set
## to get started.
setMethod("rand", "MixedTypeEngine", function(object, n, keepall = FALSE, ...) {
  ## 'chop' applies pre-existing cutpoints and factor labels to
  ## appropriate data feature-rows
  chop <- function(mtengine, dataset) {
    OC <- mtengine@cutpoints
    for (I in 1:length(OC)) {
      if (is.null(OC[[I]]$breaks | "continuous" == OC[[I]]$Type)) {
        next # do nothing for continuous data
      }
      ## othwerise, discretize following the rules
      dataset[,I] <-  cut(dataset[,I],
                          breaks = OC[[I]]$breaks,
                          labels = OC[[I]]$labels)
    }
    dataset
  }

  eng <- as(object, "CancerEngine")
  dataset <- rand(eng, n)
  hazy <- blur(object@noise, dataset)
  binned <- chop(object, hazy)
  if (keepall) {
    binned <- list(raw = dataset, noisy = hazy, binned = binned)
  }
  binned
})


## Default settings for low levels of noise in clinical data. Users can
## if they wish, go back to the original gene-expression based noise model.
ClinicalNoiseModel <- function(shape = 1.02, scale = 0.05/shape) {
  new("NoiseModel",
      additiveOffset = 0,
      additiveScale = rgamma(1, shape=shape, scale=scale),
      multiplicativeScale = 0)
}

## UNFINISHED
## nExtraBlocks = how many features do you want that are unrelated to clusters?
ClinicalBHP <- function(nFeatures,
                        nExtraBlocks = NULL,
                        meanBlockSize = NULL,
                        sigmaBlockSize =  0,
                        minBlockSize =  NULL,
                        mu0 = 6,
                        sigma0 = 1.5,
                        rate = 28.11,
                        shape = 44.25,
                        p.cor = 0.6,
                        wt.cor = 5) {
  ## TEMP: FIX THIS!
  if (is.null(nExtraBlocks)) {
    nExtraBlocks <-  3
    meanBlockSize <- 5
    minBlockSize <- 5
  }
  BlockHyperParameters(nExtraBlocks,
                       meanBlockSize,
                       sigmaBlockSize,
                       minBlockSize,
                       mu0,
                       sigma0,
                       rate,
                       shape,
                       p.cor,
                       wt.cor
                       )
}

## Original version from CEC master's thesis.
## Arguments:
##    nFeatures = number of features
##    nClusters = number of clusters/subtypes
##    weighted = logical value used to define prevalence of subtypes
##    bHyp = object of class BlockHyperParameters, if not NULL
ClinicalEngine <- function(nFeatures, nClusters, isWeighted, bHyp = NULL){
  hitfn <- ifelse(nFeatures < 15, # small feature size
                  function(n) 2,
                  ifelse(nFeatures < 45,
                         function(n) 3, # medium feature size
                         function(n) 5 # large feature size
                         )
                  )
  ## nPossible = number of possible 'hits' based on multi-hit model
  ## of cancer, where smaller numbers are well separated
  nPos <- function(nF) {
    ifelse(nF < 12,
           nF,   # no correlation; each feature will be a seaprate hit
           ifelse (nF < 50,
                   round(nF/3),
                   ifelse(nF < 100,
                          round(nF/5),
                          20)
                   )
           )
  }

  mod <- CancerModel(name="Clinical Simulation Model (Raw)",
                     nPossible = nPos(nFeatures), # number of possible hits
                     nPattern = nClusters,        # number of subtypes
                     HIT = hitfn,
                     prevalence = Prevalence(isWeighted, nClusters)
                     )
  if (is.null(bHyp)) {
    bHyp <- ClinicalBHP(nFeatures) # use defaults based on features space
  }
  eng <- makeBlockStructure(mod, bHyp)
}


### Below here, code is taken from Simulations_More.Rmd, which was used to
### crete simulation's for CEC's master's thesis.
Prevalence <- function(weighted, k){
  if(weighted=="equal"){
    a = 100
    prev <- rdirichlet(1, rep(a, k))
    as.numeric(prev)
  } else if(weighted=="unequal"){
    if(k<=8){
      a = 10
      prev <- rdirichlet(1, rep(a, k))
      as.numeric(prev)
    } else if (k>=16){
      my.min=-1
      w = 5
      while(my.min<0.01){
      alphas <- c(rep(w*1,k/4), rep(w*2,k/4), rep(w*4,k/4), rep(w*8,k/4))
      prev <- rdirichlet(1, c(alphas))
      my.min <- min(prev)
    } # end while loop
      as.numeric(prev)
    } # end if/else cluster size for unequal weights
  } # end if/else statement
} # end function

EngineBuilder <- function(p, f, k, weighted){
  if(f < 15){
    hitfn = function(n) 2
  } else {
    hitfn = function(n) 5
  }
  mod <- CancerModel(name="Cluster Simulation Model",
                     nPossible=f/3, #number of possible 'hits' based on multi-hit theory of cancer, where smaller numbers are well separated
                     nPattern=k, #number of subtypes
                     HIT = hitfn,
                     prevalence = Prevalence(weighted, k)
                     )
  bHyp <- BlockHyperParameters(nExtraBlocks=0, #how many features do you want that are not related to the clusters? 
                               meanBlockSize=3,
                               sigmaBlockSize= 0,
                               minBlockSize= 3,
                               mu0=6,
                               sigma0=1.5,
                               rate=28.11,
                               shape=44.25,
                               p.cor=0.6,
                               wt.cor=5
                               )
  eng <- makeBlockStructure(mod, bHyp)
}


## KRC: Not needed, since "MakeData <- rand" is equivalent to this definition
MakeData <- function(eng, p){
  sim <- rand(eng, p)
}

## KRC: Probably not needed; only point is that the choice of parameters
## is taken away from the user.
Noisy <- function(sim){
  mod <- ClinicalNoiseModel(30, 40, 0.10)
  noi <- blur(mod, sim)
}


# Determine if the feature is nominal or ordinal, and store that information in a vector.
catType <- sample(c(1,2), 1, replace=TRUE) #1 is ordinal, 2 is nominal



## KRC: Begining version of "chooseDataTypes", though the temptation
## to call it "Alligator" is hard to resist.
Allocator <- function(cont, bina, cat, nFeats){ #percentages or fractions of each data type
  numCont <- nFeats*cont
  numBin <- nFeats*bina
  numCat <- nFeats*cat
  feats <- sample(1:nFeats, nFeats, replace=FALSE)
  contFeats <- feats[1:numCont]
  binFeats <- feats[(numCont+1):(numCont+numBin)]
  catFeats <- feats[(numCont+numBin+1):(numCont+numBin+numCat)]
  x <- cbind(contFeats, binFeats, catFeats)
  x <- as.data.frame(x)
}#close function


